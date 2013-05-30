#!/usr/bin/env python
#
from collections import defaultdict
import sys
#
import drdcommon
from drdvcf import Vcf, VcfSnp

usage = """
tool = snp_by_func_cons.py

Given an annotated vcf file report how many snps we have
per each func consequence.

min: minimum number of samples to consider before dropping
a snp because all the genotypes are the same

Usage:
  $ cat vcf | tool <min> > stats.tsv
"""

def process_snps(vcf, min_num_samples):
  skipped, total = 0, 0
  hc = defaultdict(lambda: 0) # fc -> count
  for l in vcf.each_snp():
    snp = VcfSnp(l)
    if snp.annotated == False:
      raise(Exception('Found a snp that is not annotated: %s' % l))
    if len(snp.gtypes()) >= min_num_samples and snp.all_gtypes_the_same():
      skipped += 1
    else:
      hc[snp.func_cons] += 1
      total += 1
  hc["SAME_GENOTYPE_ON_ALL_SAMPLES_SKIPPED"] = skipped
  hc["TOTAL"] = total
  return hc

def report(h_fc_counts):
  print "func_cons\tcounts"
  for fc in h_fc_counts:
    print "%s\t%s" % (fc, h_fc_counts[fc])

def do_work(fd_vcf, min_num_samples):
  vcf = Vcf(fd_vcf)
  vcf.load_meta_header()
  report(process_snps(vcf, min_num_samples))

def main():
  if len(sys.argv) != 2:
    drdcommon.error("Wrong # of args", usage)
  if drdcommon.data_in_stdin() == False:
    drdcommon.error("Need data in stdin.", usage)

  min_num_samples = int(sys.argv[1])
  fd_vcf = drdcommon.xopen("-")
  do_work(fd_vcf, min_num_samples)
  fd_vcf.close()

if __name__ == "__main__":
  main()
