#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
import pysam
import logging
#
import drdcommon

drdcommon.setup_logging(logging)

tool = os.path.basename(__file__)
usage = """
tool = %s

Given a bam file (mrfast) and a egenotype probe file:
  . load probes in mem
  . per each alignment (including alternative XA:Z ...)
    - hits_probe?
      * probe_id +1
  . each probe
    - Id, n_hits

Usage:
  $ %s <bam> <probes.tsv>
""" % (tool, tool)

MIN_MAPQ = 10
def good_alignment(ar):
  return not ar.is_unmapped and ar.mapq >= MIN_MAPQ and not ar.is_duplicate

def gen_chrm_lenghts(samfile):
  """ Given a bam (samfile), give me list of chrms and the sizes of chrms"""
  d_genome = {}
  chrm_names = [e['SN'] for e in samfile.header['SQ']]
  for name, size in zip(chrm_names, samfile.lengths):
    d_genome[name] = size
  return chrm_names, d_genome

class Probes(object):
  def __init__(self):
    self.ids  = []      # [id1, id2....]
    self.info = {}      # [id_index] -> [chrm, coor]
    self.locations = {} # [chrm][coordiante] -> id_index
    self.n_locations = 0

  def load(self, probes_fn):
    id_index=0
    for l in open(probes_fn):
      chrm, coor, _id, five, three = l.split("\t")[0:5]
      self.ids.append(_id)
      self.info[id_index] = [chrm, coor]
      self.set_locations_for(chrm, coor, id_index, five, three)
      id_index+=1

  def set_locations_for(self, chrm, coor, id_index, five, three):
    if chrm not in self.locations:
      self.locations[chrm] = {}
    self.locations[chrm][coor] = id_index
    self.n_locations+=1
    i_coor = int(coor)
    _len_five, _len_three = len(five), len(three)
    for i in range(1, _len_five+1):
      self.locations[chrm][i_coor+i] = id_index
      self.n_locations+=1
    for i in range(1, _len_three+1):
      self.locations[chrm][i_coor-i] = id_index
      self.n_locations+=1

  def size(self):
    return len(self.info), self.n_locations

  def aln_hits_probe(self, chrm, coordinates):
    """Tell me if the alignment (chrm + coordinates) contains/hits any probe
    ------------------[XXXXXXXXXXXXXXX]-------------
                   ****************              !HIT
                         ********************    !HIT
                     ********************        HIT
    """
    hit_start, hit_end = False, False
    pindex = ""
    for c in coordinates:
      if chrm in self.locations:
        if not hit_start and c in self.locations[chrm]:
          hit_start = True
          pindex    = self.locations[chrm][c]
        if hit_start and c not in self.locations[chrm]:
          hit_end = True
    if hit_start and hit_end:
      return pindex
    return 0

  def get_id(self, pindex):
    return self.ids[pindex]

def iterate_reads(samfile, probes):
  total = 0
  hits = []
  for alignedread in samfile.fetch():
    total += 1
    t_name = samfile.getrname(alignedread.tid)
    pid = probes.aln_hits_probe(t_name, alignedread.positions)
    if pid > 0:
      hits.append(pid)
    if total % 100000 == 0:
      logging.info("%d|%d" % (total, len(hits)))
  logging.info("%d|%d" % (total, len(hits)))
  return hits

def check_input():
  if len(sys.argv) != 3:
    drdcommon.error("Wrong # of args", usage)
  bam_fn, probes_fn = sys.argv[1], sys.argv[2]
  logging.info("bam: %s probes: %s" % (bam_fn, probes_fn))
  if not os.path.isfile(bam_fn):
    drdcommon.error("Invalid bam file.", usage)
  return bam_fn, probes_fn

def dump_results(hits, probes):
  d = defaultdict(lambda: 0)
  for pindex in hits:
    d[probes.get_id(pindex)]+=1
  for pid, n_hits in d.items():
    print "%s\t %d" % (pid, n_hits)

def main():
  bam_fn, probes_fn = check_input()
  samfile = pysam.Samfile(bam_fn, "r", header=None)
  probes = Probes()
  logging.info("Loading probes ...")
  probes.load(probes_fn)
  logging.info("# of locations/loci: " + str(probes.size()))
  hits = iterate_reads(samfile, probes)
  dump_results(hits, probes)
  samfile.close()

if __name__ == "__main__":
  main()
