#!/usr/bin/env python
#
# Given the maf dataset in stdin plots a barplot of the allele freqs.
#
import matplotlib
matplotlib.use('Agg') #
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
#
import drdcommon
from drdvcf import Vcf, VcfSnp
import drdplots
import sys
from drdmath import log_it

def main():
  if len(sys.argv) == 3:
    df = pd.read_table(sys.argv[1])
    title = "MAF CRV"
    labels = ["0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5" ]
    drdplots.barplot(df.counts, labels, title, ofn=sys.argv[2])
  else:
    drdcommon.error("Wrong number of args. Need input tsv file and output png.")

if __name__ == "__main__":
  main()
