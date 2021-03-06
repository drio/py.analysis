#!/usr/bin/env python

#import drdplots
#import numpy as np
from pandas import DataFrame, read_table
#import matplotlib.pyplot as plt
#from pylab import *
from collections import defaultdict
from drdcommon import xopen
import argparse
import sys

def two_dec_points(v):
  return "{0:.2f}".format(v)

def misassembly(v, _max):
  # v: value call in all samples ; max: maximun number of # of samples with the same call
  return len(v[v==1]) >= _max or len(v[v==-1]) >= _max

# Test with:
# $ rm -f foo.txt; ../downstream.py --input_fn joined.test.csv -o foo.txt -n 300  -t 100  &&  cat foo.txt | sed  's/\t/ /g' | sort -k2,2n > f; mv f ./foo.txt
class StateMachine:
  def __init__(self, threshold):
    self.h_events        = {}
    self.event_type      = 0 # no del or ins
    self.start           = 0 # start coor of current event
    self.size            = 0
    self.th              = threshold

  def enter_new_event(self, coor, diff_rd):
    # print "  enter_new_event(%s, %s)" % (coor, diff_rd)
    self.start      = coor
    self.event_type = self.diff_to_etype(diff_rd)
    self.size       = 0
    self.prev_coor  = coor

  def diff_to_etype(self, diff_rd):
    if diff_rd >= self.th:
      return 1
    elif diff_rd <= -self.th:
      return -1
    else:
      return 0

  def save(self, coor):
    assert self.event_type != 0
    self.size += coor - self.prev_coor
    # print "  SAVE S:" , self.start, "Z:", self.size
    self.h_events[self.start] = (self.event_type, self.size)

  """
  def done(self):
    if self.event_type != 0 and self.size != 0:
      self.save()
  """

  def get_results(self):
    l_coor, l_sizes, l_types = [], [], []

    for coor, t in self.h_events.iteritems():
      l_coor.append(coor)
      l_sizes.append(t[1])
      l_types.append(t[0])

    return DataFrame({ "coor": l_coor, "size": l_sizes, "type": l_types, })

  def change_state(self, coor, diff_rd):
    if self.event_type == 0: # We are not in an event
      if self.diff_to_etype(diff_rd) != 0:
        self.enter_new_event(coor, diff_rd)
    else: # in event
      if self.diff_to_etype(diff_rd) == self.event_type:
        self.size      += coor - self.prev_coor
        self.prev_coor = coor
      else: # different type of event
        self.save(coor)
        if self.diff_to_etype(diff_rd) != 0:
          self.enter_new_event(coor, diff_rd)
        else: # no event
          self.event_type = 0

def call_stretches(df, threshold, what=1):
  sm = StateMachine(threshold)
  for coor, diff_rd in df.itertuples():
    # print "LINE: ", coor, diff_rd
    sm.change_state(coor, diff_rd)
  return sm.get_results()

def count_num_of_events(df):
  samples = defaultdict(lambda: {'i':0, 'd':0})

  for name, bins in df.iterrows():
    v = bins.values # all the RD values in all the bins for sample name
    samples[name]['i'] = len(v[v==1])
    samples[name]['d'] = len(v[v==-1])

  return DataFrame(samples).transpose()

def map_it(df, th=1, no_call=False):
  def logic(v):
    if v > th:
      return 1
    elif v < -th:
      return -1
    else:
      return 0
  return df.applymap(logic)

def filter_bins_where_all_samples_have_same_call(df, _max):
  n = 0
  for coor in df.columns:
    if misassembly(df[coor], _max):
      del df[coor]
  return df

def normalize(df, rd_ref):
  # sample_name -> [rd_values]
  h = defaultdict(lambda: [])

  for sample_name in df.columns:
    avg = np.average(df[sample_name])
    a_rd = df[sample_name]
    new_values = a_rd*(rd_ref/avg)
    h[sample_name] = new_values

  new_df = DataFrame(h, index=df.index.values)
  def logic(v):
    if log == 0:
      return 0.0001
    else:
      return v

  return new_df.applymap(logic)

def adjust_resolution(df, n_win=10):
  samples = defaultdict(lambda: [])

  for name, bins in df.iterrows():
    b_values = []
    indices, _tmp_i  = [], []
    for i, v in enumerate(bins.values):
      b_values.append(v)
      _tmp_i.append(i)
      if i % n_win == 0:
        samples[name].append(np.average(b_values))
        b_values = []
        indices.append(_tmp_i[0])
        _tmp_i   = []
    # Last bin
    if len(b_values) > 0:
      samples[name].append(np.average(b_values))
      indices.append(_tmp_i[0])

  return DataFrame(samples, index=indices).transpose()

def plot_hm(data, output_name):
  fig = plt.figure(figsize=(21,10))
  ax = fig.add_subplot(111)
  ax.set_xlabel('chrm')
  ax.set_ylabel('sample')
  #ax.set_axis_off()
  fig.add_axes(ax)
  cax = ax.imshow(data, aspect = 'normal', interpolation='nearest', cmap=cm.Blues)
  #cax = ax.imshow(data, aspect = 'normal', cmap=cm.Blues)
  cbar = fig.colorbar(cax, orientation='horizontal')
  #ax.pcolor(data)
  savefig(output_name)

def load_data(fn):
  # index: coordinates; columns: samples
  df = read_table(fn, sep=" ")
  coordinates = df['start']
  del df['chrm']; del df['start']
  h = {}
  for name, values in df.transpose().iterrows():
    h[name] = []
    for i, v in enumerate(values):
      h[name].append(v)

  return DataFrame(h, index=coordinates)

def compute_dist_event_freqs(d_calls):
  freq  = defaultdict(lambda : 0)
  sizes = []
  for sample, d_coors in d_calls.iteritems():
    for coor, t in d_coors.iteritems():
      if t[0] != 0:
        sizes.append(t[1])
        freq[coor] += 1

"""
  fig = plt.figure()
  ax = fig.add_subplot(111)
  v  = [i for i in freq.values() if i != 1]
  plt.xlabel('# of samples')
  plt.ylabel('Frequency')
  plt.title(r'Distribution of CNV event frequencies (del 1)')
  ax.hist(v, 50, facecolor='green')
  savefig("sample.freq.png")

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.xlabel('size')
  plt.ylabel('Frequency')
  plt.title(r'Distribution of CNV event sizes (2,100)')
  v  = [i for i in sizes if i > 2 and i < 100]
  ax.hist(v, 100, facecolor='green')
  savefig("sizes.freq.png")
"""

def parse_args():
  parser = argparse.ArgumentParser(description='cnv pipeline')
  parser.add_argument('-i', '--input_fn', metavar='input_fn', required=True,
                        dest='input_fn', action='store',
                        help='input data file')

  parser.add_argument('-o', '--output_fn', metavar='output_fn', required=False,
                        dest='output_fn', action='store',
                        help='output data file')

  parser.add_argument('-r', '--resolution', metavar='resolution', required=False,
                        dest='change resolution to this value', action='store', type=int,
                        help='Change the resolution of the resulting')

  parser.add_argument('-t', '--threshold', metavar='threshold', required=True,
                        dest='threshold', action='store', type=float,
                        help='read depth threashold for calling an event')

  args = parser.parse_args()
  args.input_fn = xopen(args.input_fn)
  if not args.output_fn:
    args.output_fn = sys.stdout
  return args

def run():
  args = parse_args()

  # chrm, start, n_reads_in_sample1, n_reads_in_sample2, ......
  df = load_data(args.input_fn)

  df = call_stretches(df, args.threshold)
  df.to_csv(args.output_fn, sep="\t")

run()


