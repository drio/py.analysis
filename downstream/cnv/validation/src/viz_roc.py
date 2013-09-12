#!/usr/bin/env python

#import drdplots
#import numpy as np
# from pandas import DataFrame, read_table
from matplotlib.pyplot import *
#from pylab import *
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
        if th not in data:
          data[th] = {}
        data[th][sen] = fdr
  return data

def do_work():
  fig       = figure(figsize=(14,12))
  ax        = fig.add_subplot(1, 1, 1)
  for th, h2 in load_data().iteritems():
      fdrs, sens = [], []
      for s in sorted(h2.iterkeys()): # get all the (fdr, sens) for a particular th
          sens.append(s)
          fdrs.append(h2[s])
      # have data(fdrs, sens) for a truth set; and particular threshold
      ax.set_title("Performance CNV caller")
      ax.set_ylabel("Sensitivity")
      ax.set_xlabel("False Discovery Rate")
      ax.plot(fdrs, sens,'o', markersize=10, label=th)
      ax.legend(loc='best', shadow=True)
  savefig('roc.png', dpi=200, bbox_inches='tight')

matplotlib.rcParams.update({'font.size': 16, 'family': 'bold'})
do_work()



