#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os, re
#
import drdcommon

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Generate windows (reads) based on illumina reads.

$ samtools view 18277.bam |\
  head -1000 | ./windify.py |\
  java -jar /stornext/snfs6/rogers/drio_scratch/bb/local/picard/FastqToSam.jar |\
  SAMPLE_NAME=XXXX QUALITY_FORMAT=Illumina FASTQ=/dev/stdin OUTPUT=out.bam

""" % (tool_name)

def dump_fq_record(_id, seq):
  print "@%s" % _id
  print seq
  print "+"
  print "J" * len(seq)

def windify(_id, seq, w_size):
  windows = {}
  i = 0
  while i <= len(seq)-w_size:
    w_seq = seq[i:i+w_size]
    w_id  = _id + "_w" + str(i)
    windows[w_id] = w_seq
    i+=1
  return windows

def do_work(fd_reads):
  for line in fd_reads:
    s = line.split()
    _id, seq, qual = s[0], s[9], s[10]
    for w_id, w_seq in windify(_id, seq, 31).items():
      dump_fq_record(w_id, w_seq)

def main():
  if len(sys.argv) == 1:
    fd_reads = drdcommon.xopen("-")
    do_work(fd_reads)
    fd_reads.close()
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
