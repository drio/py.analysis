#!/usr/bin/env python

from collections import defaultdict
import drdcommon
import argparse, sys
import random

class BitMask(object):
  def __init__(self):
    self.mask = 0

  def set(self, pos):
    self.mask = self.mask | (1 << pos)

  def exists(self, pos):
    return (self.mask & (1 << pos)) > 0

def next_coor(prev_end):
  return prev_end + random.randint(10000, 20000)

def next_size():
  return random.randint(100, 20000)

def in_n_region(start, end, ref):
  num_ns, s, e = 0, int(start), int(end)
  _max = e-s/8
  for i in range(s, e+1):
    if ref.exists(i):
      num_ns += 1
      if num_ns > _max:
        return True
  return False

if len(sys.argv) < 4 or not drdcommon.data_in_stdin():
  sys.stderr.write("cat ref.fa | tool <n_events> <chrm> <chrm_size>" + "\n")
  sys.exit(1)

_ = sys.argv
n_events, chrm, chrm_size = int(_[1]), _[2], int(_[3])

# Store N locations
sys.stderr.write("Loading N locations ..." + "\n")
ref = BitMask()
i = 0
for l in drdcommon.xopen("-"):
  if l[0] != '-':
    for c in l.rstrip():
      if c.upper() == 'N':
        ref.set(i)
    i += 1

# Generate events
sys.stderr.write("Generating events ..." + "\n")
coor = next_coor(100)
i = 0
while (i < n_events):
  # chrm start end 0..n
  s = next_size()
  if coor + s < chrm_size:
    if i % 2 == 0: # deletion
      start, end, _type = coor, coor + s, 0
    else:
      start, end, _type = coor, coor + s, str(random.randint(1, 10))
  else:
    break

  if not in_n_region(start, end, ref):
    i += 1
    print chrm, start, end, _type

  coor = next_coor(end)
