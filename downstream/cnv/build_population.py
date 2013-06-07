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
tool = %s

Look for files matching the pattern in the current directory, load them
and dump the population level cnv calls.

Usage:
  $ %s <pattern>
""" % (tool, tool)

def load_data(sid, fn, h):
  for l in open(fn):
    # chrm pos type #_bins
    chrm, coor, etype, n_bins = l.split()
    chrm = re.sub(r'(^[cC]hrm?)', '', chrm)
    if etype == 'ins' or etype == 'del': # TODO: Fix bug in previous tool
      h[etype][chrm][int(coor)][sid] = int(n_bins)

def main():
  if len(sys.argv) != 3:
    drdcommon.error("Wrong # of args", usage)
  pattern  = sys.argv[1]
  re_id    = sys.argv[2] # Regular expression to extract id
  l        = lambda:defaultdict(l)
  h        = l()
  l_ids    = []

  for fn in drdcommon.files_in_dir(".", pattern):
    try:
      sid = re.search(re_id, fn).group(1)
      l_ids.append(sid)
    except:
      raise(Exception('Problems extracting id using regular expression.'))
    load_data(sid, fn, h)

  for t, one in h.items():
    for chrm, two in one.items():
      for coor, three in two.items():
        # TODO print ids
        sys.stdout.write("%s %s %s %s " % (str(chrm), str(coor), t, len(three)))
        for sid in l_ids:
          if sid in three:
            sys.stdout.write(str(three[sid]) + " ")
          else:
            sys.stdout.write("0 ")
        print ""

if __name__ == "__main__":
  main()
