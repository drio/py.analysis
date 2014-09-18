#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
#
import drdcommon

tool = os.path.basename(__file__)
usage = """
tool = %s

Usage:
  $ %s <polyphen_calls> <shift_calls>
""" % (tool, tool)

def load_predictions(i_file, chrm_col, coor_col, columns):
    pre = {}
    header = True
    for l in drdcommon.xopen(i_file):
        if header:
            header=False
            continue
        s = l.split("\t")
        key = drdcommon.canonic_chrm(s[chrm_col]) + "_" + s[coor_col]
        _tmp = []
        for c in columns:
            _tmp.append(s[c])
        pre[key] = "\t".join(_tmp)
    return pre

def merge(pp, si, num_of_cols):
    r = {}
    for k, v in pp.iteritems():
        r[k] = {}
        r[k]["pp"] = v

    for k, v in si.iteritems():
        if k not in r:
            r[k] = {}
        r[k]["si"] = v

    for k, h in r.iteritems():
        line = []
        line.append(k.replace('_', "\t"))
        for si_or_pp in [ "pp", "si" ]:
            if si_or_pp in h:
                line.append(h[si_or_pp])
            else:
                line.append("\t".join((["-"] * num_of_cols)))
        print "\t".join(line)

def main():
    if len(sys.argv) == 3:
        pp_calls, sift_calls = sys.argv[-2:]

        int_cols = [11, 59]
        chrm_col, coor_col = -6, -5
        pre_pp = load_predictions(pp_calls, chrm_col, coor_col, int_cols)

        int_cols = [16, 25]
        chrm_col, coor_col = 2, 3
        pre_si = load_predictions(sift_calls, chrm_col, coor_col, int_cols)

        merge(pre_pp, pre_si, len(int_cols))
    else:
        print usage

if __name__ == '__main__':
    main()
