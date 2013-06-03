#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
#
import drdcommon
from common import check_sorted, good_alignment, reads_at_locus

tool = os.path.basename(__file__)
usage = """
tool = %s

Given the log2ratios, call stretches of cnv events.
We will iterate over the log2ratios and count how many
adjacent events we have that are above (or below) the
threshold.

Usage:
  $ cat lo2ratios | %s <log2ratio threshold> calls.txt
""" % (tool, tool)

class CnvStateMachine(object):
  def __init__(self, ratios_stream, threshold):
    self.threshold     = threshold
    self.ratios_stream = ratios_stream

    self.num_bins      = 0
    self.event_type    = None
    self.locus         = None

  def an_insertion(self, log_ratio):
    return log_ratio > self.threshold

  def a_deletion(self, log_ratio):
    return log_ratio < -self.threshold

  def report_strech(self):
    if self.num_bins > 0:
      print "%s %s %s %s" % (self.locus[0], self.locus[1], self.event_type , self.num_bins)

  def entering_strech(self, curr_type):
    return not self.event_type and curr_type

  def streching(self, curr_type):
    return self.event_type and curr_type == self.event_type

  def change_of_event_while_in_event(self, curr_type):
    return self.num_bins > 0 and curr_type != self.event_type

  def event_ended(self, curr_type):
    return self.event_type and not curr_type

  def change_state(self, curr_type, chrm, start):
    if self.entering_strech(curr_type):
      self.locus      = (chrm, start)
      self.num_bins   = 1
      self.event_type = curr_type
    elif self.streching(curr_type):
      self.num_bins = self.num_bins + 1
    elif self.change_of_event_while_in_event(curr_type):
      self.report_strech()
      self.locus      = (chrm, start)
      self.num_bins   = 1
      self.event_type = curr_type
    elif self.event_ended(curr_type):
      self.report_strech()
      self.num_bins   = 0
      self.event_type = None
    else: # Strech of nothing
      pass

  def run(self):
    for line in self.ratios_stream:
      # Chr1 1 5199 1000 1282 0.358
      chrm, start, end, n_reads_ref, n_reads, log_ratio = line.split()
      if self.an_insertion(float(log_ratio)):
        self.change_state("ins", chrm, start)
      elif self.a_deletion(float(log_ratio)):
        self.change_state("del", chrm, start)
      else:
        self.change_state(None, chrm, start)
    self.report_strech() # There may be something to report

def main():
  if len(sys.argv) != 2:
    drdcommon.error("Wrong # of args", usage)
  if not drdcommon.data_in_stdin():
    drdcommon.error("No data in stdin.", usage)
  ratios_stream = drdcommon.xopen("-")
  threshold = float(sys.argv[1])
  CnvStateMachine(ratios_stream, threshold).run()

if __name__ == "__main__":
  main()
