#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
import pysam
#
import drdcommon
from common import check_sorted, good_alignment, gen_chrm_lenghts, reads_at_locus

tool = os.path.basename(__file__)
usage = """
tool = %s

Computes the bins for a set of alignments.

Usage:
  $ %s <num_reads_per_window> <bam> <target_chrm>

Knowing your coverage, the desired bin size and your read length we can
find out the num_reads_per_window:

  num_reads_per_window = x_coverage * (window_size/read_length)
""" % (tool, tool)

def print_bins(samfile, max_n_reads, bam_name, chrm, chrm_name_to_length):
  start = 1
  bin_reads = set([])
  for pucol in samfile.pileup(chrm, 1, chrm_name_to_length[chrm]):
    bin_reads = reads_at_locus(pucol, bin_reads)
    if len(bin_reads) >= max_n_reads:
      end = pucol.pos
      print "%s %d %d %d" % (chrm, start, end, len(bin_reads))
      n_bin_reads = 0
      start = end
      bin_reads = set([])

def main():
  if len(sys.argv) != 4:
    drdcommon.error("Wrong # of args", usage)
  n_reads, bam_name, target_chrm = int(sys.argv[1]), sys.argv[2], sys.argv[3]
  if not os.path.isfile(bam_name):
    drdcommon.error("Invalid bam file.", usage)
  samfile = pysam.Samfile(bam_name, "rb")
  chrms, chrm_name_to_length = gen_chrm_lenghts(samfile)
  if not target_chrm in chrm_name_to_length:
    samfile.close()
    drdcommon.error("Chrm not present in bam header.", usage)
  print_bins(samfile, n_reads, bam_name, target_chrm, chrm_name_to_length)
  samfile.close()

if __name__ == "__main__":
  main()
