#!/usr/bin/env python

import sys
import os
import glob
import re
from pprint import pprint as pp
import re
import pandas as pd
import matplotlib.pyplot as plt
import drdcommon as drd
import json

def compute_dist(_files):
    df_c = pd.read_table(_files[0])[1:]
    df_o = pd.read_table(_files[1])[1:]
    df_c = df_c[(~df_c['#CHROM'].str.contains('random|chrX|chrY|chrM'))]
    df_o = df_o[(~df_o['#CHROM'].str.contains('random|chrX|chrY|chrM'))]
    o = pd.merge(df_o, df_c, on=["#CHROM", "START", "END", "GC%"])
    o = o[o['IS_CONTROL'] == 'Y']
    o.head()
    return o.COPYNUMBER.round(2).value_counts()


def plot_dist(freqs):
    plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    #axis([0,5, 0, 200000])
    plt.xlabel("copynumber")
    plt.ylabel("frequency")
    plt.title("Distribution of CNV values for the control regions")
    plot(freqs.index, freqs, "ow")
    plt.savefig('dist.cnv.control.png', format='pdf')


def extract_from_log(fn):
    d = {}
    comp = re.compile("^(\w+) Average Read Depth: ([\d\.]+), Standard Deviation: ([\d\.]+)$")
    for l in open(fn):
        m = re.match(comp, l)
        if m:
            win, avg, std = m.groups()
            d[win] = {"avg": avg, "std": std}
    return d 


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print "Usage: tool <cw_norm.bed> <copynumber.bed> <output.log> <sample_id>"
        sys.exit(1)
    else:
        _files = sys.argv[1:]
        freqs = compute_dist(_files[0:2])

        cana = extract_from_log(_files[-2])
        cana["id"] = _files[-1]
        cana["dist"] = freqs.to_dict()
        print json.dumps(cana, sort_keys=True, indent=4)
        #print freqs
