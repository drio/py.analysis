#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, re, os.path
#
import drdcommon
from common import check_sorted, good_alignment, reads_at_locus

tool = os.path.basename(__file__)
usage = """
tool = %s; Join all the ratios files.

Usage:
  $ %s <pattern> <regex_for_id> <files_dir>
""" % (tool, tool)

def load_data(sid, fn, h):
  for l in drdcommon.xopen(fn):
    # chrm start end n_reads n_reads_ref log2ratio
    chrm, start, end, nref, nr, log = l.split()
    chrm = re.sub(r'(^[cC]hrm?)', '', chrm)
    h[chrm][int(start)][sid] = int(nr)

def main():
  if len(sys.argv) != 4:
    drdcommon.error("Wrong # of args", usage)
  pattern  = sys.argv[1]
  re_id    = sys.argv[2] # Regular expression to extract id
  f_dir    = sys.argv[3]

  l_ids    = []
  l        = lambda:defaultdict(l)
  h        = l() # hold all data in mem

  for fn in drdcommon.files_in_dir(f_dir, pattern):
    try:
      sid = re.search(re_id, fn).group(1)
      l_ids.append(sid)
    except:
      raise(Exception('Problems extracting id using regular expression.'))
    load_data(sid, fn, h)

  # print header
  sys.stdout.write("chrm start ")
  for _id in l_ids:
    sys.stdout.write("%s " % _id)
  print ""

  for chrm, one in h.items():
    for start, two in one.items():
      sys.stdout.write("%s %s " % (str(chrm), str(start)))
      for sid, nr in two.items():
        sys.stdout.write(str(nr) + " ")
      print ""

if __name__ == "__main__":
  main()
