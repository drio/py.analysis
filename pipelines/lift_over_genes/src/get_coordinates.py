#!/usr/bin/env python

import drdcommon
import sys
import re

def help(msg):
  sys.stderr.write("ERROR: " + msg + "\n")
  sys.stderr.write("Usage: cat gtf.txt | tool list_genes_names.txt > genes.coor.bed\n")
  sys.exit(1)

# Main
if not drdcommon.data_in_stdin():
  help("Need data in stdin")

if len(sys.argv) != 2:
  help("Invalid list of arguments")

gene_names = {}
for l in drdcommon.xopen(sys.argv[1]):
  name = l.split()[0]
  gene_names[name] = True
drdcommon.log("%s genes loaded." % len(gene_names))

found = 0
for l in drdcommon.xopen("-"):
  s = l.split("\t")
  if s[2] == "CDS":
    chrm, start, end, _list = s[0], s[3], s[4], s[8]
    name = None
    for e in _list.split(";"):
      _ = e.split()
      if len(_) == 2 and _[0] == "gene_name":
        name = re.sub('\"', '', _[1])
        if name in gene_names:
          print "chr%s %s %s %s" % (chrm, start, end, name)
          found += 1
          del gene_names[name]

drdcommon.log("%s genes found." % found)
