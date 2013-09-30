#!/usr/bin/env python

import drdcommon
import sys

def help(msg):
  sys.stderr.write("ERROR: " + msg + "\n")
  sys.stderr.write("Usage: cat refgene.txt | tool list_genes_names.txt > genes.coor.bed\n")
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
  s = l.split()
  chrm, start, end, name = s[2], s[4], s[5], s[12]
  if name in gene_names:
    print "%s %s %s %s" % (chrm, start, end, name)
    found += 1

drdcommon.log("%s genes found." % found)
