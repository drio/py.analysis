#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os, re, fnmatch
#
import drdcommon
from drdvcf import Vcf, VcfSnp
from drdmath import std_dev

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Given probe alignments and eg hits tell me if the probe
is specific enough.

$ %s <pattern for files>
""" % (tool_name, tool_name)

def do_work(pattern):
  data = []
  for f in drdcommon.files_in_dir('.', pattern):
    fd = open(f)
    for l in fd.readline():
      print l
    fd.close()
def main():
  if len(sys.argv) == 2:
    pattern = sys.argv[1]
    do_work(pattern)
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
