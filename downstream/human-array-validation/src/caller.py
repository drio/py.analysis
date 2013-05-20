#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os, re
#
import drdcommon
from drdvcf import Vcf, VcfSnp
from drdmath import std_dev

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Given probe alignments and eg hits tell me if the probe
is specific enough.

$ %s <sam_stream> <eg_hits_stream> <min_mapq> <pm_hits> > calls.txt

pm_hits: +/- value to use when determine if eg hits passes filter.
""" % (tool_name, tool_name)

class ValidateChip(object):
  def __init__(self, fn_sam, fn_hits, min_mq, pm_hits):
    self.fn_sam  = fn_sam
    self.fn_hits = fn_hits
    self.min_mq  = min_mq
    self.probe_info = {}
    self.pm_hits = pm_hits

  def do_work(self):
    self.process_alignments()
    most_common_n_hits, max_times = self.find_most_common_cov_val()
    self.process_hits(most_common_n_hits)
    self.make_calls()

  # HWI-ST115:330:C0HRPACXX:2:2315:16354:36604      177     Chr1    6       37      101M    ChrU    153732487       0       GGAGGTTGAGACCAG
  def process_alignments(self):
    drdcommon.log("Processing alignments")
    init_doesnt_pass_eg_hits = False
    for line in drdcommon.xopen(self.fn_sam):
      if line[0] != '@':
        s = line.split()
        probe_id, mq = s[0], int(s[4])
        has_good_qual = mq > self.min_mq
        self.probe_info[probe_id] = [ has_good_qual, init_doesnt_pass_eg_hits ]

  def find_most_common_cov_val(self):
    drdcommon.log("Finding common coverage value.")
    d = defaultdict(int)
    for n_hits, p_id in self.iterate_over_eg_cov():
      d[n_hits] += 1
    n_hits_most_freq, max_times = 0, 0
    for n_hits, times in d.items():
      if times > max_times:
        max_times, most_comm_n_hits = times, n_hits
    drdcommon.log("n_hits:%s times:%s" % (most_comm_n_hits, max_times))
    return most_comm_n_hits, max_times

  def process_hits(self, most_common_n_hits):
    drdcommon.log("Processing hits")
    _min, _max = most_common_n_hits - self.pm_hits, most_common_n_hits + self.pm_hits
    for n_hits, p_id in self.iterate_over_eg_cov():
      has_resonable_coverage = n_hits >= _min and n_hits <= _max
      self.probe_info[p_id][1] = has_resonable_coverage

  def make_calls(self):
    for _id, a in self.probe_info.items():
      ok_qual, ok_hits = a[0], a[1]
      call = "P"
      if not ok_qual and not ok_hits:
        call = "B"
      elif not ok_qual:
        call = "M"
      elif not ok_hits:
        call = "E"
      print "%s\t%s" % (_id, call)

  def iterate_over_eg_cov(self):
    fd_hits = drdcommon.xopen(self.fn_hits)
    for line in fd_hits:
      s = line.split()
      n_hits, p_id = int(s[0]), s[1].rstrip()
      yield n_hits, p_id
    fd_hits.close()

def main():
  if len(sys.argv) == 5:
    fn_sam, fn_hits   = sys.argv[1:3]
    min_mapq, pm_hits = [ int(i) for i in sys.argv[3:]]
    ValidateChip(fn_sam, fn_hits, min_mapq, pm_hits).do_work()
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
