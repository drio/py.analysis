#!/usr/bin/env python
#
import matplotlib
matplotlib.use('Agg')
import pandas as pd
from collections import defaultdict
#
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
#
import drdcommon
from drdvcf import Vcf, VcfSnp
import drdplots
from drdmath import log_it

def plot(ofn, x, y, title, xlabel, ylabel):
  #x_ticks = [i for i in range(1,int(max(x))+5) if i % 5000 == 0 or i == 1]
  #y_ticks = [-3, -2, -1, 0, 1, 2, 3]
  dot_size=1

  fig = plt.figure()
  fig.set_size_inches(20,3)
  ax = fig.add_subplot(1, 1, 1)
  #ax.add_line(Line2D([0.5, 0.5], [0, 1], transform=ax.transAxes, linewidth=2, color='b'))
  l = ax.axhline(y=0, color="red", linewidth=1)
  #ax.set_xticks(ticks=x_ticks)
  #ax.set_xticklabels(x_ticks)
  #ax.set_yticks(ticks=y_ticks)
  #ax.set_yticklabels(y_ticks)
  ax.scatter(x, y, dot_size)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  #plt.grid(True)
  plt.title(title)
  plt.savefig(ofn, dpi=100, bbox_inches='tight')

def process_data(fd):
  """ chrm start end control_n_reads sample_n_reads log2ratio
      Chr1 1 5642 1500 997 -0.589
  """
  logs = []
  for line in fd:
    ratio = float(line.rstrip().split()[-1])
    #clear_event = ratio > 0.2 or ratio < -0.2
    #if clear_event:
    logs.append(ratio)
    #else:
    #  logs.append(0)
  return logs

def main():
  if len(sys.argv) == 3:
    logratios = process_data(drdcommon.xopen("-"))
    bin_nums  = range(1, len(logratios)+1)
    title     = sys.argv[1]
    output_fn = sys.argv[2]
    plot(output_fn,
      bin_nums, logratios, title, xlabel="bin #", ylabel="log2ratios (sample/control)")
  else:
    drdcommon.error("Wrong number of args. <title> <output.filename>")

if __name__ == "__main__":
  main()
