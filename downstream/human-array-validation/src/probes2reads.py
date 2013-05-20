#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os, re
#
import drdcommon
from drdvcf import Vcf, VcfSnp
from drdmath import std_dev

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Given a probe file (stdin), generate a fasta file.
""" % (tool_name)

VALID_SEQ = re.compile(r"^[AGCTN]+$")

def valid_seq(seq):
  match = VALID_SEQ.search(seq)
  if match:
    return True
  else:
    return False

# 13      77361753        rs8000613       AATAGCATTGATAGG AATCGACCTAAGTGT A       G       0.38825214899713467     0.465616045845272206    0.146131805157593123
# TODO: drop probes with non nucleotide values
def do_work(fd_probes):
  for line in fd_probes:
    s = line.split()
    _id = s[2]
    seq = s[3] + 'N' + s[4]
    if valid_seq(seq):
      print "@%s" % _id
      print "%s" % seq
      print "+"
      print "J" * len(seq)
    else:
      sys.stderr.write("Skipping %s\n")

def main():
  if len(sys.argv) == 1:
    fd_probes = drdcommon.xopen("-")
    do_work(fd_probes)
    fd_probes.close()
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
