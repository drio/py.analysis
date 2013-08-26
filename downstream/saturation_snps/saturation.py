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
  def __init__(self, stream):
    self.stream    = stream
    self.n_samples = 1
    self.subs      = defaultdict(lambda : 0)
    self.indels    = defaultdict(lambda : 0)
    self.counts    = {}
    self.more_samples_to_process = True
    self.__doWork()

  def csv(self, sep=","):
    csv = "num_samples" + sep + "substitutions" + sep + "indels" + "\n"
    index_subs, index_indels = 0, 1
    for n_samples, counts in self.counts.items():
      csv += "%s%s%s%s%s\n" % (n_samples, sep, str(counts[index_subs]), sep, str(counts[index_indels]))
    return csv

  def __doWork(self):
    while self.more_samples_to_process:
      self.__jump_header()
      self.__process_snps()
      self.__update_counts()
      self.n_samples += 1
    return self.counts

  def __update_counts(self):
    n_subs   = self.__seen_in_more_than_n_samples(self.subs, AT_LEAST_SEEN)
    n_indels = self.__seen_in_more_than_n_samples(self.indels, AT_LEAST_SEEN)
    #self.counts[self.n_samples] = (len(self.subs), len(self.indels))
    self.counts[self.n_samples] = (n_subs, n_indels)
    sys.stderr.write(">> #:%s SUBS:%s INDELS:%s MEM(Mbytes):%s\n" % \
      (self.n_samples, n_subs, n_indels, drdcommon.memory_usage()))

  def __seen_in_more_than_n_samples(self, d, n):
    num_snps = 0
    for coor, n_samples_with_snp in d.items():
      if n_samples_with_snp >= n:
        num_snps += 1
    return num_snps

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
    self.more_samples_to_process = False

  def __in_header(self, l):
    return l[0] == '#'

def parse_args():
  parser = argparse.ArgumentParser(description='staturation capture stats')

  parser.add_argument('-g', '--genes', metavar='genes', required=False,
                        dest='calls_chrm', action='store', default=False,
                        help='Report genes discovered instead of snps.')

  args = parser.parse_args()
  return args

def main():
  args = parse_args()
  stream = drdcommon.xopen("-")
  if not drdcommon.data_in_stdin():
    drdcommon.error(usage)
  print Saturation(stream).csv("\t")
  stream.close()

if __name__ == "__main__":
  main()
