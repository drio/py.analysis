#!/usr/bin/env python

"""
Given a bed file (stdin) and a min/max average read depth
parameter, filter genes where any of the exons does not
match the range of RD mean.
"""

from collections import defaultdict
import drdcommon
import argparse, sys

if len(sys.argv) != 3:
  sys.stderr.write("cat my.bed | tool <min> <max> > in 2> out" + "\n")
  sys.exit(1)

out = sys.stdout.write
err = sys.stderr.write

_ = sys.argv
a_min, a_max = float(_[1]), float(_[2])

#chrm    start   end     gene    exon_number     transcript_number       count   mean    std     min     25%     50%     75%     max
#Chr7    21807454        21807531        76P     1       76P-201 78.0    116.47  14.53   0.0     116.0   119.5   122.0   127.0

# Per each gene, save the average RD of each of its exons
genes = {}
first_line = True
for l in drdcommon.xopen("-"):
  if first_line:
    first_line = False
    continue

  s = l.strip().split()
  chrm, start, end, g_name = s[0:4]
  start, end = int(start), int(end)
  mean = float(s[7])

  if g_name not in genes:
    genes[g_name] = { "chrm" : chrm,
                      "start": start,
                      "end"  : end,
                      "means": [] }

  # Update start/end as necessary
  cg = genes[g_name] # current gene
  if start < cg["start"]:
    cg["start"] = start
  if end > cg["end"]:
    cg["end"] = end

  cg["means"].append(mean)

# Filter and dump
for name, v in genes.items():
  dump = out
  for m in v["means"]:
    if m < a_min or m > a_max:
      dump = err
      break
  dump("%s\t%s\t%s\t%s\n" % (v["chrm"], v["start"], v["end"], name) )
