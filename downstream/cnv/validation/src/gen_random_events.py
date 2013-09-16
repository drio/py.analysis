#!/usr/bin/env python

from collections import defaultdict
import drdcommon
import argparse, sys
import random

if len(sys.argv) < 4:
  sys.stderr.write("tool <n_events> <chrm> <chrm_size>" + "\n")
  sys.exit(1)

def next_coor(prev_end):
  return prev_end + random.randint(100, 1000)

def next_size():
  return random.randint(100, 20000)

_ = sys.argv
n_events, chrm, chrm_size = int(_[1]), _[2], int(_[3])

coor = next_coor(100)
for i in range(0, n_events):
  # chrm start end 0..n
  s = next_size()
  if coor + s < chrm_size:
    if i % 2 == 0: # deletion
      start, end, _type = coor, coor + s, 0
    else:
      start, end, _type = coor, coor + s, str(random.randint(1, 10))
  else:
    break

  coor = next_coor(end)
  print chrm, start, end, _type

