#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os.path
import drdcommon
import argparse, sys

def parse_args():
  parser = argparse.ArgumentParser(description='cnv pipeline')
  parser.add_argument('-e', '--events', metavar='events', required=True,
                        dest='events', action='store',
                        help='List of events to introduce in the genome')

  parser.add_argument('-r', '--reference', metavar='reference', required=True,
                        dest='reference', action='store',
                        help='original fasta file of the reference')

  args = parser.parse_args()
  args.events_stream    = drdcommon.xopen(args.events)
  args.reference_stream = drdcommon.xopen(args.reference)
  return args

def done(args):
  args.events_stream.close()
  args.reference_stream.close()

class BitMask(object):
  def __init__(self):
    self.mask = defaultdict(lambda: 0)

  def set(self, chrm, pos):
    self.mask[chrm] = self.mask[chrm] | (1 << pos)

  def exists(self, chrm, pos):
    return (self.mask[chrm] & (1 << pos)) > 0

def load_events(stream):
  bmask = BitMask()
  count = 0
  events = defaultdict(lambda: {})
  for l in stream:
    if not l[0] == '#':
      sl = l.split()
      chrm, coor = sl[0], int(sl[1])
      is_deletion  = len(sl) == 4 and sl[3] == '0'
      is_insertion = len(sl) == 5 and sl[1] == sl[2]

      if not is_deletion and not is_insertion:
        raise Exception('Incorrect entry event: ' + l)
      if chrm in events and coor in events[chrm]:
        raise Exception('Duplicated event')
      if bmask.exists(chrm, coor):
        raise Exception('Event already there.')

      if is_deletion:
        events[chrm][coor] = [ int(sl[1]), int(sl[2]) ] # start, end
      if is_insertion:
        events[chrm][coor] = [ int(sl[1]), int(sl[3]), int(sl[4]) ] # start, size, times
      bmask.set(chrm, coor)
      count += 1

  sys.stderr.write(">> %d events loaded\n" % count)

def run():
  args = parse_args()
  bed = load_events(args.events_stream)
  done(args)

run()

