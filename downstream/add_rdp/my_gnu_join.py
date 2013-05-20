#!/usr/bin/env python
#
import logging, os
from drdcommon import *
from drdvcf import Vcf, VcfSnp

tool = os.path.basename(__file__)
logging.basicConfig(format=tool + ': %(asctime)s >> %(message)s', level=logging.DEBUG)

def read_new_lines(fds):
  lines = [fd.readline().rstrip() for fd in fds] # read one line from all the files
  lines = [l for l in lines if l != '']          # deal with EOF
  return lines

def process():
  fds = [] # file descriptors
  for fn in sys.argv[1:]:
    fds.append(xopen(fn))

  delim = "\t"
  lines = read_new_lines(fds)
  while lines:
    chrm, coor, cov = lines[0].split()
    sys.stdout.write(chrm + "_" + coor + delim + cov + delim) # we add chrm coor and coverage for first sample
    for l in lines[1:]:
      sys.stdout.write(l.split()[2] + delim)
    print
    lines = read_new_lines(fds)

  for fd in fds:
    fd.close()

def main():
  if len(sys.argv) > 1:
    process()
  else:
    error("usage: my_gnu_join file1 file2 ... fileN")

if __name__ == "__main__":
  main()
