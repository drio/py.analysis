#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os.path, argparse
import pysam
#
import drdcommon
from drdvcf import VcfSnp

tool = os.path.basename(__file__)
usage = """
tool = %s

Given a stream of vcf files (stdin) computes the increment of snps as a function
of the samples processed.

-g: compute how many new genes we are discovering when more samples are added.

Usage:
  $ gzip -cd *.vcf.gz | %s > staturation.txt

NOTE: We are only counting events seen in at least 2 animals.

""" % (tool, tool)

MIN_QUAL = 20
AT_LEAST_SEEN = 2

class Saturation(object):
    def __init__(self, stream, at_least_seen):
        self.stream    = stream
        self.n_samples = 1
        self.at_least_seen = at_least_seen
        self.subs      = defaultdict(lambda : 0)
        self.indels    = defaultdict(lambda : 0)
        self.genes     = defaultdict(lambda : 0)
        self.genes_partial = {}
        self.counts    = {}
        self.more_samples_to_process = True
        self.__doWork()

    def csv(self, sep=","):
        csv = "num_samples" + sep + "substitutions" + sep + "indels" + sep + "genes" + "\n"
        index_subs, index_indels, index_genes = 0, 1, 2
        for n_samples, counts in self.counts.items():
            csv += "%s%s%s%s%s%s%s\n" % (n_samples, sep,
                                     str(counts[index_subs]), sep,
                                     str(counts[index_indels]), sep,
                                     str(counts[index_genes]))
        return csv

    def __doWork(self):
        while self.more_samples_to_process:
            self.__jump_header()
            self.__process_snps()
            self.__update_counts()
            self.n_samples += 1
        return self.counts

    def __update_counts(self):
        n_subs   = self.__seen_in_more_than_n_samples(self.subs, self.at_least_seen)
        n_indels = self.__seen_in_more_than_n_samples(self.indels, self.at_least_seen)
        n_genes  = self.__seen_in_more_than_n_samples_for_genes(self.genes_partial, self.at_least_seen)
        self.counts[self.n_samples] = (n_subs, n_indels, n_genes)
        sys.stderr.write(">> #:%s SUBS:%s INDELS:%s GENES:%s MEM(Mbytes):%s AT_LEAST_SEEN:%s\n" % \
          (self.n_samples, n_subs, n_indels, n_genes, drdcommon.memory_usage(), self.at_least_seen))

    def __seen_in_more_than_n_samples(self, d, n):
        num_snps = 0
        for coor, n_samples_with_snp in d.items():
            if n_samples_with_snp >= n:
                num_snps += 1
        return num_snps

    def __seen_in_more_than_n_samples_for_genes(self, d, n):
        for g_name, _ in d.items():
            self.genes[g_name] += 1
        self.genes_partial = {} # Start empty for the next sample.
        return len([gn for gn, count in self.genes.items() if count >= n])

    def __jump_header(self):
        for l in self.stream:
            if self.__is_header_line(l):
                break

    def __is_header_line(self, l):
        return l[0] == '#' and l[1] != '#'

    def __process_snps(self):
        for l in self.stream:
            if self.__in_header(l):
                self.more_samples_to_process = True
                return
            else:
                snp = VcfSnp(l)
                if snp.has_high_quality(MIN_QUAL):
                    if snp.is_a_substitution():
                        self.subs[snp.coordinate()] += 1
                    else:
                        self.indels[snp.coordinate()] += 1
                    if snp.annotated:
                        self.genes_partial[snp.gene] = True

        self.more_samples_to_process = False

    def __in_header(self, l):
        return l[0] == '#'


def parse_args():
    parser = argparse.ArgumentParser(description='staturation capture stats')

    parser.add_argument('-g', '--genes', metavar='genes', required=False,
                          dest='genes', action='store', default=False,
                          help='Report genes discovered instead of snps.')

    parser.add_argument('-a', '--at_least_seen', metavar='at_least_seen', required=False,
                          dest='at_least_seen', action='store', default=2, type=int,
                          help='Count snp when seen in at least that number of samples.')


    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    stream = drdcommon.xopen("-")
    if not drdcommon.data_in_stdin():
        drdcommon.error(usage)
    print Saturation(stream, args.at_least_seen).csv("\t")
    stream.close()

if __name__ == "__main__":
    main()
