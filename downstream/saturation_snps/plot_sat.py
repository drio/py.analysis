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
import matplotlib.pyplot as plt
from drdmath import log_it

_usage = """
Usage:
  cat ../saturation/old/all/staturation.all.txt  | awk '{print $1" "$2}' | grep -v num |\
  plot_stat.py "title" "x_label" "y_label" > output.png
"""

def plot(x, y, title, xlabel, ylabel):
    #x_ticks = [i for i in range(1,int(max(x))+5) if i % 5 == 0 or i == 1]
    dot_size=12

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #ax.set_xticks(ticks=x_ticks)
    ax.scatter(x, y, dot_size)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.title(title)
    plt.savefig(sys.stdout, bbox_inches='tight')

def process_data(fd):
    x, y = [], []
    for line in fd:
        s = line.rstrip().split()
        if len(s) == 2:
            _id, count = s
            x.append(_id)
            y.append(count)
    return x, y

def main():
    if not drdcommon.data_in_stdin():
        drdcommon.error("I need a data stream in stdin.", usage=_usage)
    if not len(sys.argv) == 4:
        drdcommon.error("Wrong number of parameters", usage=_usage)

    title, _xl, _yl = sys.argv[1:]
    x, y = process_data(drdcommon.xopen("-"))
    plot(x, y, title, xlabel=_xl, ylabel=_yl)

if __name__ == "__main__":
    main()
