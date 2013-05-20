#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
import pysam
import math
#
import drdcommon
from common import check_sorted, good_alignment, reads_at_locus

tool = os.path.basename(__file__)
usage = """
tool = %s

Computes the log2 ratios of the number of reads seen per
each of the control window.

Usage:
  $ cat windows.txt | %s <bam>
""" % (tool, tool)

def log_2_ratio(w_n_reads, sample_n_reads):
  if sample_n_reads == 0:
    sample_n_reads = 0.00001
  _log = math.log(sample_n_reads/w_n_reads, 2)
  return str(round(_log, 3))

def compute_ratios(windows, bam_name):
  samfile = pysam.Samfile(bam_name, "rb")
  for l in windows:
    w_chrm, w_start, w_end, w_n_reads = l.split()
    reads = set([])
    for pucol in samfile.pileup(w_chrm, int(w_start), int(w_end)+1):
      reads = reads_at_locus(pucol, reads)
    print l.rstrip() + " " + str(len(reads)) + " " + log_2_ratio(int(w_n_reads), len(reads))
  samfile.close()

def main():
  if len(sys.argv) != 2:
    drdcommon.error("Wrong # of args", usage)
  if not drdcommon.data_in_stdin():
    drdcommon.error("No data in stdin.", usage)
  windows = drdcommon.xopen("-")
  bam_name = sys.argv[1]
  if not os.path.isfile(bam_name):
    drdcommon.error("Invalid bam file.", usage)
  compute_ratios(windows, bam_name)

if __name__ == "__main__":
  main()
