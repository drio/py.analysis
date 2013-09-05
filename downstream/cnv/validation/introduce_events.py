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
      chrm = sl[0]
      start, end, times = [int(v) for idx, v in enumerate(sl) if idx > 0]

      is_deletion  = len(sl) == 4 and sl[3] == '0'
      is_insertion = len(sl) == 4 and sl[3] != '0'
      if not is_deletion and not is_insertion:
        raise Exception('Incorrect entry event: ' + l)
      if chrm in events and start in events[chrm]:
        raise Exception('Duplicated event')
      if bmask.exists(chrm, start):
        raise Exception('Event already there.')

      events[chrm][start] = [ end, times ]
      bmask.set(chrm, start)
      count += 1

  sys.stderr.write(">> %d events loaded\n" % count)
  return events

def gen_new_reference(events, ref_stream):
  first_chrm = True
  in_del = 0
  ins_chunk, ins_size, ins_times = [], 0, 0
  c_chrm = ""
  o_pt   = -1 # pointer in the original reference

  for line in ref_stream:
    if line[0] == '>':
      c_chrm = line[1:-1]
      o_pt = 0
      if not first_chrm:
        print
      print line.rstrip()
      first_chrm = False
      continue

    for c in line:
      if c == "\n": continue

      if o_pt != 0 and o_pt % 20 == 0: print

      o_pt += 1
      if ins_size == 0 and in_del == 0: # we are not dealing with an event
        if o_pt in events[c_chrm]: # we have to introduce an event here
          entering_deletion = events[c_chrm][o_pt][1] == 0
          if entering_deletion:
            in_del = events[c_chrm][o_pt][0] - o_pt

          entering_insertion = events[c_chrm][o_pt][1] > 0
          if entering_insertion:
            ins_size, ins_times =  (events[c_chrm][o_pt][0] - o_pt) + 1, events[c_chrm][o_pt][1] + 1

      #
      if in_del > 0:
        in_del -= 1
        continue

      if ins_size > 0:
        ins_chunk.append(c)
        ins_size -= 1
        if ins_size == 0:
          sys.stdout.write((''.join(ins_chunk) * ins_times ))
          ins_chunk = []
        continue

      sys.stdout.write(c)
  print

def run():
  args = parse_args()
  bed = load_events(args.events_stream)
  gen_new_reference(bed, args.reference_stream)
  done(args)

run()

