#!/usr/bin/env python
#
import sys, os
from drdcommon import *

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Given a bed output from samtools depth and a original locations.bed,
add the missing locations from the samtools depth output. Samtools does
not report a location if the bam does not have any coverage.

Usage:
  $ samtools depth -b locations.bed my.bam | %s locations.bed > fixed_depth.bed
""" % (tool_name, tool_name)

def do_work(fdd, fdl): # fd depth, fd locations
  ldep, lloc = fdd.readline(), fdl.readline()

  while lloc:
    if ldep:
      dep_chrm, dep_coor = ldep.split()[0:2]
    loc_chrm, loc_start, loc_end = lloc.split()

    if dep_chrm == loc_chrm and dep_coor == loc_end:
      sys.stdout.write(ldep)
      ldep = fdd.readline()
    else:
      sys.stdout.write("%s\t%s\t%s\n" % (loc_chrm, loc_end, 0))
    lloc = fdl.readline()

def main():
  if len(sys.argv) == 2 and data_in_stdin():
    fd_depth     = xopen("-")
    fd_locations = xopen(sys.argv[1])
    do_work(fd_depth, fd_locations)
    fd_depth.close()
    fd_locations.close()
  else:
    if not data_in_stdin():
      error("No data in stdin.", usage)
    else:
      error("Wrong # of params", usage)

if __name__ == "__main__":
  main()

