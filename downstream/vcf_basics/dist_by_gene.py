#!/usr/bin/env python
#
from __future__ import division
import pandas as pd
from collections import defaultdict
import sys
#
import drdcommon
from drdvcf import Vcf, VcfSnp
import drdplots as drdplots
from drdmath import std_dev

usage = """
tool = dist_by_func_cons

Given an annotated vcf file report how many snps we have
per each func consequence.

Usage:
  $ cat vcf | tool 2> dist.tsv
"""
def process_snps(vcf):
  h_genes = defaultdict(lambda: 0)
  for l in vcf.each_snp():
    snp = VcfSnp(l)
    if snp.annotated == False:
      raise(Exception('Found a snp that is not annotated: %s' % l))
    if snp.impact == "HIGH" and snp.gene != "":
      h_genes[snp.gene] += 1
  return h_genes

def report(h_genes):
  print "gene_name\tcounts"
  for g in h_genes:
    print "%s\t%s" % (g, h_genes[g])

def do_work(fd_vcf):
  vcf = Vcf(fd_vcf)
  report(process_snps(vcf))

def main():
  if len(sys.argv) == 2:
    drdcommon.error("Wrong # of args", usage)
  if drdcommon.data_in_stdin() == False:
    drdcommon.error("Need data in stdin.", usage)

  fd_vcf = drdcommon.xopen("-")
  do_work(fd_vcf)
  fd_vcf.close()


if __name__ == "__main__":
  main()
