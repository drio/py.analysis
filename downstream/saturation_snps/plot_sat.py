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

def plot(ofn, x, y, title, xlabel, ylabel):
  x_ticks = [i for i in range(1,int(max(x))+5) if i % 5 == 0 or i == 1]
  dot_size=12

  fig = plt.figure()
  ax = fig.add_subplot(1, 1, 1)
  ax.set_xticks(ticks=x_ticks)
  ax.scatter(x, y, dot_size)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.grid(True)
  plt.title(title)
  plt.savefig(ofn, dpi=400, bbox_inches='tight')

def process_data(fd):
  h = {'n_samples': [], 'n_subs': [], 'n_indels': [] }
  first_line = True
  for line in fd:
    if not first_line and line != "\n":
      n_samples, n_subs, n_indels  = line.rstrip().split('\t')
      h['n_samples'].append(float(n_samples))
      h['n_subs'].append(float(n_subs))
      h['n_indels'].append(float(n_indels))
    else:
      first_line = False
  return pd.DataFrame(h)

def main():
  if len(sys.argv) == 1:
    fd = drdcommon.xopen("-")
    df = process_data(fd)
    title="Saturation of substitutions CRV WGS"
    plot("crv.subs.saturation.png", \
      df['n_samples'], df['n_subs'], title, xlabel="# samples", ylabel="# uniq substitutions")

    title="Saturation of indels CRV WGS"
    plot("crv.indels.saturation.png", \
      df['n_samples'], df['n_indels'], title, xlabel="# samples", ylabel="# uniq indels")

    fd.close()
  else:
    drdcommon.error("Wrong number of args. Just need std values in stdin.")

if __name__ == "__main__":
  main()
