#!/usr/bin/env python
#
# scatter plot var allele frequencies
#
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *

def doit(fn):
  data = pd.read_csv(fn,names=['counts', 'naa'], sep='\s+')
  total = 0.
  for e in data.counts.values: total+=e
  f = lambda x: x*100/total
  data['percentage'] = data.counts.apply(f)
  return data

wgs = doit('alt.wgs.txt')
wes = doit('alt.wes.txt')
