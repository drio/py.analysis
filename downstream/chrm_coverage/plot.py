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

def plot(x, y, title, xlabel, ylabel):
    #x_ticks = [i for i in range(1,int(max(x))+5) if i % 5 == 0 or i == 1]
    #dot_size=12
    fig = plt.figure(figsize=(16,6))
    ax = fig.add_subplot(1, 1, 1)
    #ax.set_xticks(ticks=x_ticks)
    #ax.scatter(x, y, dot_size)
    ax.scatter(x, y)
    ax.set_ylim([0,100])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.title(title)
    plt.savefig(sys.stdout, bbox_inches='tight')

def process_data(fd, window_size=10000):
    x, y, win_vals = [], [], []
    i = 0
    w_num = 1
    for line in fd:
        s = line.rstrip().split()
        if len(s) == 3:
            chrm, coor, value = s
            win_vals.append(int(value))
            i += 1
            if i == 10000:
                x.append(w_num)
                y.append(sum(win_vals) / len(win_vals))
                w_num += 1
                i = 0
                win_vals = []

    if len(win_vals) > 0:
        x.append(sum(win_vals) / len(win_vals))
        y.append(w_num)

    return x, y


def main():
    if not drdcommon.data_in_stdin():
        drdcommon.error("I need a data stream in stdin.", usage="-")
    if not len(sys.argv) == 2:
        drdcommon.error("Wrong number of parameters", usage="-")

    title = sys.argv[1]
    x, y = process_data(drdcommon.xopen("-"))
    plot(x, y, title, xlabel="genomic window", ylabel="Average Read Depth")

if __name__ == "__main__":
    main()
