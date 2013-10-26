#!/usr/bin/env python

"""
Given a set of targets/exons and the bases coverage for them, report
generate another bed with the desc stats of each target/exon.
"""

from collections import defaultdict
import drdcommon
import argparse, sys
import random
from pandas import Series

if len(sys.argv) != 3:
  sys.stderr.write("tool <target/exons bed file> <bed base coverage>" + "\n")
  sys.exit(1)

_ = sys.argv
fn_targets, fn_base_cov = _[1], _[2]

# Load exons (targets)
depth = {}
for t in drdcommon.xopen(fn_targets):
    l = t.strip()
    sl = [c for c in l.split()] # splitted line
    chrm, start, end = sl[0:3]
    if chrm not in depth:
        depth[chrm] = {}
    if start in depth[chrm]:
        raise(Exception('Two exons starting in same location! bailing out: ' + l))
    for i in range(int(start), int(end)+1):
        depth[chrm][i] = 0

# Fill Read depth per base
for t in drdcommon.xopen(fn_base_cov):
    chrm, coor, rd = t.strip().split()
    depth[chrm][int(coor)] = int(rd)

out = sys.stdout.write

# header
header = ["chrm", "start", "end", "gene", "exon_number", "transcript_number"]
h_stats = list(Series([1,2,3]).describe().index)
print "\t".join(header + h_stats)

# Compute stats
for t in drdcommon.xopen(fn_targets):
    # DRY
    l = t.strip()
    sl = [c for c in l.split()] # splitted line
    chrm, start, end = sl[0:3]

    out("\t".join(sl))

    l = []
    for i in range(int(start), int(end)+1):
         l.append(int(depth[chrm][i]))

    for i in Series(l).describe():
        out("\t" + str(round(i, 1)))

    print

