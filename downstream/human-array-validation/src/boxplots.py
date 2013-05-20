#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os, re, fnmatch
#
import drdcommon
from drdvcf import Vcf, VcfSnp
from drdmath import std_dev

import matplotlib
matplotlib.use('Agg') # hack to avoid DISPLAY error when in the console
import matplotlib.pyplot as plt

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Generate a box plot for all the distribution values in
files that match the pattern file.

$ %s <pattern for files>
""" % (tool_name, tool_name)

def boxplot(data, title="title here", ofn="boxplot.png", y_limit=None):
  fig = plt.figure()
  fig.set_size_inches(14,4)
  ax = fig.add_subplot(1,1,1)
  ax.set_xticklabels(range(len(data)), rotation=45)
  if y_limit:
    plt.ylim(0, y_limit)
  ax.set_title(title)
  ax.boxplot(data)
  plt.savefig(ofn, dpi=150)

def load_data(pattern):
  data = []
  for f in drdcommon.files_in_dir('.', pattern):
    sys.stderr.write(f + "\n")
    l_cov_vals = []
    for l in open(f):
      l_cov_vals.append(int(l.split()[0]))
    data.append(l_cov_vals)
  if len(data) == 1:
    data = data[0]
  sys.stderr.write(str(len(data)) + "\n")
  return data

def main():
  if len(sys.argv) == 2:
    pattern = sys.argv[1]
    data = load_data(pattern)
    boxplot(data, title=pattern, y_limit=50)

  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
