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

for l in drdcommon.xopen("-"):
  s = l.split("\t")
  if s[2] == "CDS":
    chrm, start, end, _list = s[0], s[3], s[4], s[8]
    g_name, e_name, t_name = None, None, None

    for e in _list.split(";"):
      _ = e.split()
      if len(_) == 2 and _[0] == "transcript_name":
        t_name = re.sub('\"', '', _[1])
      if len(_) == 2 and _[0] == "gene_name":
        g_name = re.sub('\"', '', _[1])
      if len(_) == 2 and _[0] == "exon_number":
        e_num = re.sub('\"', '', _[1])
      if g_name in gene_names and g_name and e_num and t_name:
        print "chr%s %s %s %s %s %s" % (chrm, start, end, g_name, e_num, t_name)
        break
