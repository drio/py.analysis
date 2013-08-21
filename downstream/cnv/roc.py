#!/usr/bin/env python

#import drdplots
#import numpy as np
# from pandas import DataFrame, read_table
#import matplotlib.pyplot as plt
#from pylab import *
from collections import defaultdict
import drdcommon
import argparse, sys

def process_calls(ds_t, calls_f, min_event_size):
  fp, tp, fn, sensitivity, fdr = 0, 0, 0, 0, 0

  for c in calls_f:
    line = c.split()
    if len(line) < 3:
      raise Exception('Expecting: chrm, start, end, ...', c)
    chrm, start, end = line[0:3]
    chrm, start, end = drdcommon.canonic_chrm(chrm), int(start), int(end)
    if end-start >= min_event_size:
      ds_t.query(chrm, start, end)

  tp, fp, fn = ds_t.compute_metrics(chrm)
  sys.stderr.write(">> tp: %s fp: %s fn: %s" % (tp, fp, fn) + "\n")
  # of the total true events, how many we find with our tool
  sensitivity = float(tp) / (tp + fn)
  # of the total our system call, how many are not in the true event set
  fdr         = float(fp) / (fp + tp)

  return sensitivity, fdr

class Truth(object):
  """
  Data structure to hold the truth (hits from other tools)
  dict(chrm) = list of coordinates where each cell is the number of tools that
               called that event
  """
  def __init__(self, buffer_size): # Use low size for testing
    self.ds            = defaultdict(lambda: {})
    self.event_hits    = []

    self.buff_size     = buffer_size
    self.fp            = 0 # false positives
    self.tp            = 0 # true positives
    self.n_true_events = 0

  def load(self, truth_f, min_event_size, calls_chrm):
    for l in truth_f:
      line = l.split()
      if len(line) != 4:
        raise Exception("Expecting lines with: chrm start size id. Got: ", l)

      chrm, start, end = self.process_line(line)

      if end-start > min_event_size:
        if calls_chrm:
          if calls_chrm == chrm:
            self.add_event(chrm, start, end)
        else:
          self.add_event(chrm, start, end)

    sys.stderr.write(">> total # of true events loaded: %s\n" % self.n_true_events)
    return self

  def process_line(self, line):
    chrm, start, end, _ = line
    return drdcommon.canonic_chrm(chrm), int(start), int(end)

  def add_event(self, chrm, start, end):
    self.n_true_events += 1
    self.event_hits.append(0)
    index_new_entry = len(self.event_hits)-1
    for c in xrange(start-self.buff_size, end+1+self.buff_size, 1):
      if c > 0:
        self.ds[chrm][c] = index_new_entry

  def query(self, chrm, start, end):
    """ Is there any event in the truth that is being hit? """
    for coor in range(start, end+1):
      coor_hits_event = chrm in self.ds and coor in self.ds[chrm]
      if coor_hits_event:
        event_index = self.ds[chrm][coor]
        if self.event_hits[event_index] == 0: # No previous hits
          self.event_hits[event_index] = 1
          self.tp += 1
        return True

    self.fp += 1
    return False

  def compute_metrics(self, chrm):
    return self.tp, self.fp, self.n_true_events - self.tp

def test_ds_t(ds_t):
  assert ds_t.query('1', 200) == 2
  assert ds_t.query('1', 201) == 2
  assert ds_t.query('1', 199) == 2
  assert ds_t.query('1', 204) == 1
  assert ds_t.query('1', 196) == 1
  assert ds_t.query('1', 500) == 0

def parse_args():
  parser = argparse.ArgumentParser(description='cnv pipeline')
  parser.add_argument('-c', '--calls', metavar='calls', required=True,
                        dest='calls_f', action='store',
                        help='input bed file with the calls of your classifier')

  parser.add_argument('-t', '--truth', metavar='truth', required=True,
                        dest='truth_f', action='store',
                        help='list of validated cnv calls')

  parser.add_argument('-o', '--output_f', metavar='output_fn', required=False,
                        dest='output_f', action='store', default=sys.stdout,
                        help='list of validated cnv calls')

  parser.add_argument('-w', '--buffer', metavar='win_buffer', required=False,
                        dest='buffer_size', action='store', type=int, default=0,
                        help='buffer to use when checking hit against truth')

  parser.add_argument('-m', '--min_size', metavar='min_size', required=False,
                        dest='min_size', action='store', type=int, default=0,
                        help='min size of the events you want to consider when loading the truth. ')

  parser.add_argument('-r', '--chrm', metavar='calls_chrm', required=False,
                        dest='calls_chrm', action='store',
                        help='Chrm we used for the calls.')

  args = parser.parse_args()
  args.calls_f = drdcommon.xopen(args.calls_f)
  args.truth_f = drdcommon.xopen(args.truth_f)
  return args

def run():
  args = parse_args()
  print args
  sys.stderr.write(">> Loading the truth" + "\n")
  ds_t = Truth(args.buffer_size).load(args.truth_f, args.min_size, args.calls_chrm)

  sys.stderr.write(">> Computing metrics" + "\n")
  sensitivity, fdr = process_calls(ds_t, args.calls_f, args.min_size)
  print "sensitivity fdr"
  print sensitivity, fdr

  """
  # ../roc.py -c /etc/passwd -t small.bed  -w 1
  test_ds_t(ds_t)
  """

run()


