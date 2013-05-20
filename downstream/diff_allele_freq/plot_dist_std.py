#!/usr/bin/env python
#
# std dev of var allele ratios for all the groups
#
import pandas as pd
from collections import defaultdict
#
import drdcommon
from drdvcf import Vcf, VcfSnp
import drdplots
import sys
from drdmath import log_it

def process_data(fd):
  std, counts = [], []
  for line in fd:
    c, v = line.split()
    std.append(float(v))
    counts.append(float(c))
  return std, counts

def main():
  if len(sys.argv) == 1:
    fd = drdcommon.xopen("-")
    std, counts = process_data(fd)
    title = "std dev freq of var allele ratios"

    drdplots.scatter_plot("std.dist.png",
                          std, log_it(counts, 10),
                          title=title, xlabel="std deviation",
                          ylabel="log10(counts)", dot_size=10)
    fd.close()
  else:
    drdcommon.error("Wrong number of args. Just need std values in stdin.")

if __name__ == "__main__":
  main()
