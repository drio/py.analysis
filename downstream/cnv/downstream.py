#!/usr/bin/env python

import drdplots
import numpy as np
from pandas import *
import matplotlib.pyplot as plt
from pylab import *
from collections import defaultdict

def misassembly(v, _max):
  # v: value call in all samples ; max: maximun number of # of samples with the same call
  return len(v[v==1]) >= _max or len(v[v==-1]) >= _max

class StateMachine:
  def __init__(self):
    self.h_events   = {}
    self.event_type = 0 # no del or ins
    self.start      = 0 # start coor of current event
    self.size       = 0

  def enter_new_event(self, coor, new_e_type):
    assert new_e_type != 0
    self.start      = coor
    self.event_type = new_e_type
    self.size       = 1

  def save(self):
    assert self.event_type != 0
    self.h_events[self.start] = (self.event_type, self.size)

  def done(self):
    if self.event_type != 0 and self.size != 0:
      self.save()

  def change_state(self, coor, new_e_type):
    if self.event_type == 0:
      if new_e_type != 0:
        self.enter_new_event(coor, new_e_type)
    else: # in event
      if new_e_type == self.event_type: # a strech of same event
        self.size += 1
      else: # different type of event
        self.save()
        if new_e_type != 0:
          self.enter_new_event()
        else: # no event
          self.event_type = 0

def call_it(df):
  """
  df:     a dataframe where the indices are the samples, the columns are the bin coordinates and the
          cells contain a event type
  output: dict[sample][coor] = (event_type, size)
  """
  h_calls = {}
  for sample_name, bins in df.iterrows():
    sm = StateMachine()
    [ sm.change_state(coor, val) for coor, val in zip(df.columns, bins.values)]
    sm.done()
    h_calls[sample_name] = sm.h_events
  return h_calls

def count_num_of_events(df):
  samples = defaultdict(lambda: {'i':0, 'd':0})

  for name, bins in df.iterrows():
    v = bins.values # all the RD values in all the bins for sample name
    samples[name]['i'] = len(v[v==1])
    samples[name]['d'] = len(v[v==-1])

  return DataFrame(samples).transpose()

def map_it(df, th=1):
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
  samples = defaultdict(lambda: [])

  for name, bins in df.iterrows():
    avg = np.average(bins.values)
    new_values = bins.values*(rd_ref/avg)
    new_values = new_values.astype(int)
    samples[name] = new_values

  return DataFrame(samples).transpose()

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

  return DataFrame(h, index=coordinates).transpose()

def compute_dist_event_freqs(d_calls):
  freq  = defaultdict(lambda : 0)
  sizes = []
  for sample, d_coors in d_calls.iteritems():
    for coor, t in d_coors.iteritems():
      if t[0] != 0:
        sizes.append(t[1])
        freq[coor] += 1

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

def run(fn):
  df = load_data(fn)
  df = normalize(df, 300)
  df = adjust_resolution(df, 10)
  df = np.log2(df/300)
  df = map_it(df, th=2)

  """
  print "# bins after adjusting resolution: %d " % len(df.columns)
  df = filter_bins_where_all_samples_have_same_call(df, 50)
  print "# bins after dropping bins where all samples have same call: %d" % len(df.columns)
  plot_hm(df, 'headmap.png')
  df.to_csv("calls.txt", sep="\t")
  """

  d_calls = call_it(df)
  compute_dist_event_freqs(d_calls)

  """
  df = count_num_of_events(df)
  df = df.drop('35084')
  drdplots.barplot(df['i'], df['i'].index, 'insertions', 'ins.png', x_inches=20, y_inches=5, rota=90)
  drdplots.barplot(df['d'], df['d'].index, 'deletions', 'del.png', x_inches=20, y_inches=5, rota=90)
  df.to_csv("counts.txt", sep="\t")
  """

run("joined.35087.1k.300.depth.csv")
#run("small.join.csv")
