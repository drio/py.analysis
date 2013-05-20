#!/usr/bin/env python
#
# Extract variant allele freq from input (vcf records)
#
import sys, re

RE = re.compile(r"DP4=(\d+),(\d+),(\d+),(\d+);")

def match(line):
  match = RE.search(line)
  if match:
    ref_count, var_count = 0., 0.
    ref_count += int(match.group(1))
    ref_count += int(match.group(2))
    var_count += int(match.group(3))
    var_count += int(match.group(4))
    return str(round(var_count/(ref_count+var_count), 2))
  else:
    raise("Couldn't match DP4")

def main():
  for line in sys.stdin:
    if line[0] == '#': continue
    print match(line)

if __name__=='__main__':
  sys.exit(main())

