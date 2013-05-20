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

Join calls from different samples (and probe files)

$ %s <sam_file> <patterm_eg_hits> > join.affy.tsv

$ for s in affy6 Omni exon; do %s *$s*.sam "*counts*$s*.txt" > join.$s.tsv ;done

""" % (tool_name, tool_name, tool_name)


class JoinData(object):
  ID = re.compile(r"^(\d+)\.")  # match id

  def __init__(self, sam_fn, pattern_eg_hits):
    self.pattern = pattern_eg_hits
    self.sam_fn  = sam_fn
    self.d_mq    = {}
    self.d_hits  = defaultdict(lambda: {})
    self.ids     = [] # order of samples in output file

  def do_work(self):
    self.load_mapq()
    self.load_ids_order()
    self.load_hits()
    self.dump()

  def load_mapq(self):
    drdcommon.log("Loading mapq")
    for line in drdcommon.xopen(self.sam_fn):
      if line[0] != '@':
        s = line.split()
        probe_id, chrm, coor, mq = s[0], s[2], s[3], s[4]
        self.d_mq[probe_id] = [ chrm, coor, mq ]

  def load_ids_order(self):
    for f in drdcommon.files_in_dir('.', self.pattern):
     self.ids.append(self.extract_id(f))
    drdcommon.log("# of samples loaded: %s" % len(self.ids))

  def load_hits(self):
    drdcommon.log("Loading hits")
    for f in drdcommon.files_in_dir('.', self.pattern):
      sample_id = self.extract_id(f)
      drdcommon.log("fn: %s | id: %s" % (f, sample_id))
      for line in drdcommon.xopen(f):
        s = line.split()
        assert len(s) == 2
        n_hits, p_id = s[0], s[1].rstrip()
        self.d_hits[p_id][sample_id] = n_hits

  def dump(self):
    header = ["probe_id", "chrm", "pos" , "mapq" ] + self.ids
    print "\t".join(header)

    for p_id, l_mq in self.d_mq.items(): # probes
      line = [ p_id ] + l_mq # first part of the output line
      for sample_id in self.ids: # hits in samples
        we_have_hits_on_probe_for_sample = self.d_hits[p_id].has_key(sample_id)
        n_hits = self.d_hits[p_id][sample_id] if we_have_hits_on_probe_for_sample else "0"
        line.append(n_hits)
      print "\t".join(line)

  def extract_id(self, file_name):
    match = self.ID.search(file_name)
    if match:
      return match.group(1)
    else:
      raise Exception('Could not find id in file name: %s' % file_name)

def main():
  if len(sys.argv) == 3:
    sam_fn, pattern_fn_hits = sys.argv[1:]
    JoinData(sam_fn, pattern_fn_hits).do_work()
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
