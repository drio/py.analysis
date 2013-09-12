#!/usr/bin/env python

#import drdplots
#import numpy as np
# from pandas import DataFrame, read_table
from matplotlib.pyplot import *
from pylab import *
from collections import defaultdict
import drdcommon
import argparse, sys, glob

def iter_lines(re_file, with_fn=False):
  splitted_lines = []
  for f in glob.glob(re_file):
    for l in open(f):
      s = l.split()
      _cell = (f, s) if with_fn else s
      splitted_lines.append(_cell)
  return splitted_lines

def load_data():
  data = {}
  for s in iter_lines('roc.data*.txt'):
    if len(s) == 4:
      min, th, sen, fdr = s # TODO: we are not using the min var
      data[float(fdr)] = float(sen)
  return data

def do_work():
  fig  = figure(figsize=(14,12))
  ax   = fig.add_subplot(1, 1, 1)
  fdrs, sens = [], []
  data = load_data()
  for f in sorted(data.keys()):
    fdrs.append(f)
    sens.append(data[f])
  print fdrs, sens

  ax.plot(fdrs, sens, 'o-', markersize=5)
  ax.set_title("Performance CNV caller")
  ax.set_ylabel("Sensitivity")
  ax.set_xlabel("False Discovery Rate")
  xlim([0,1])
  ylim([0,1])
  #for t, f, s in zip(thds, fdrs, sens):
  #  ax.annotate(t, xy=(f, s))
  savefig('roc.png', dpi=250, bbox_inches='tight')

matplotlib.rcParams.update({'font.size': 16, 'family': 'bold'})
do_work()



