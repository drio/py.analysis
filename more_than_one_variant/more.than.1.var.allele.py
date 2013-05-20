#!/usr/bin/env python
#
# Of the subs that have more than 1 variant allele, what's the
# ratio of the other alleles compare to all the allele seen?
#
#
#Chr1    39861   .       T       G,A     148.73  .       AC1=1;AC=2,1;AF1=0.5;AN=6;DP4=28,40,31,26;DP=160;EFF=INTERGENIC(MODIFIER|||||||||);FQ=17.1;MQ=50;PV4=0.098,1,0.0001,0.0027;SF=3,5,7;VDB=0.0392;EFF=INTERGENIC(MODIFIER|||||||||)
#GT:GQ:PL        .       .       .       0/2:47:44,.,.,0,.,255   .       0/1:99:255,0,242,.,.,.  .       0/1:99:237,0,220,.,.,.  .       .
#
import sys
import pandas as pd
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
from pylab import show
from drdcommon import *

def process_line(line):
  if line[0] == '#':
    return

  start_col_gts, alt_all_col = 9, 4
  s = line.split()
  h = {0:0, 1:0, 2:0, 3:0}
  for gt in s[start_col_gts:]:
    if gt != '.':
      try:
        for allele in [int(gt[0]), int(gt[2])]:
          h[allele]+=1
      except:
        error("Error indexing the var array with line: %s" % (line))
  return h

def process_data(fn='-'):
  fd = xopen(fn)

  total = []
  data = {0:[], 1:[], 2:[], 3:[]}
  for line in fd:
    h = process_line(line)
    t = 0.
    if h != None:
      for k,v in h.items():
        data[k].append(v)
        t += v
      total.append(t)

  if fd != '-': fd.close()
  df = DataFrame(data)
  df['total'] = total
  return df

def only(d, num_aa): # data,  num of alternative alleles
  reduced = []
  for e in d:
    if len(e) == num_aa: reduced.append(e)
  return reduced

def comp_ratio(df):
  df['var_ratio'] = (df[2] + df[3]) / df['total']
  return df

def plot_it(df, i, ofn, title="", xlabel="", ylabel=""):
  plt.subplots(1)

  plt.scatter(df.index, df['var_ratio'], 1)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.grid(True)
  plt.title(title)

  #fig.subplots_adjust(wspace=1, hspace=0.3)
  if in_ipython():
    plt.figure(i)
    print "Use show(%d) to display plot. " % i
  else:
    plt.savefig(ofn, dpi=400, bbox_inches='tight')

#files = [ 'more.wes.vcf.gz', 'more.wgs.vcf.gz']
files = [ 'small.vcf']
plots = []
for i, fn in enumerate(files):
  df = comp_ratio(process_data(fn))
  plots.append(plot_it(df, i, fn + ".png", fn + ": Alternative allele ratio (1 alt allele)", "index", "ratio var/total"))

