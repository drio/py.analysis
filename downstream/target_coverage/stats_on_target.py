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

out = sys.stdout.write
err = sys.stderr.write

_ = sys.argv
fn_targets, fn_base_cov = _[1], _[2]

err("Loading exons/targets\n")
depth = {}
n = 0
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
    n += 1
err("%s\n" % n)

err("Reading read depth bed\n")
total = n
hits = 0
for t in drdcommon.xopen(fn_base_cov):
    chrm, coor, rd = t.strip().split()
    coor = int(coor)
    total += 1
    if chrm in depth and coor in depth[chrm]:
        rd = int(rd)
        depth[chrm][coor] = rd
        if rd > 0:
            hits += 1
    if total % 500000 == 0:
        err("hits: %s total: %s \n" % (hits, total) )

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
        out("\t" + str(round(i, 2)))

    print

