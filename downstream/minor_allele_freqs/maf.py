#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
import pysam
#
import drdcommon
from drdvcf import Vcf, VcfSnp

tool = os.path.basename(__file__)
usage = """
tool = %s

Given a population vcf file (stdin) computes the minor allele frequencies.

Usage:
  $ gzip -cd merged.vcf.gz | %s > maf.tsv

""" % (tool, tool)

MIN_QUAL = 20

def process_snps(vcf):
  """ per each snp, compute the ratio of the min allele freq / total number
      of alleles seen.
      Ignore cases were we are seeing more than two alleles.
      Find how many chrmosomes you observe with the minor allele
      compute the ratio of that and the total number of alleles seen.
  """
  vcf.load_meta_header()
  total_n_chrms = vcf.num_of_samples * 2
  mafs = defaultdict(lambda: 0)
  for l in vcf.each_snp():
    snp = VcfSnp(l)
    if snp.has_high_quality(MIN_QUAL):
      a_counts = snp.alternative_allele_counts()
      a_total  = snp.total_num_alleles()
      if len(a_counts) == 1: # only 1 alternative allele
        n_chrms_with_alt_allele = a_counts[0]
        n_chrms_with_ref_allele = a_total - n_chrms_with_alt_allele
        if n_chrms_with_ref_allele <= n_chrms_with_alt_allele:
          mafs[round(n_chrms_with_ref_allele/total_n_chrms, 2)] += 1
        else:
          mafs[round(n_chrms_with_alt_allele/total_n_chrms, 2)] += 1
  return mafs

def report(mafs):
  labels = [round(i*0.05,2) for i in xrange(0,11,1)]
  labels.remove(0)
  d = defaultdict(lambda: 0)
  for m, c in mafs.items():
    for l in labels:
      if m <= l:
        d[l]+=c
        break
  print "maf\tcounts"
  mafs = d.keys()
  mafs.sort()
  for m in mafs:
    print "%s\t%s" % (m, d[m])

def do_work(fd_vcf):
  vcf = Vcf(fd_vcf)
  report(process_snps(vcf))

def main():
  if len(sys.argv) != 1:
    drdcommon.error("Wrong # of args", usage)
  if not drdcommon.data_in_stdin():
    drdcommon.error("No data in stdin.", usage)
  fd_vcf = drdcommon.xopen("-")
  do_work(fd_vcf)
  fd_vcf.close()

if __name__ == "__main__":
  main()
