#!/usr/bin/env python
#
# show how many snps we have per each func cons type
#
import sys, re

MIN_QUAL = 200
MAX_COV  = 0
MIN_COV  = 0
SF_RE    = re.compile(r"SF=([\d,]+);")
COV_RE   = re.compile(r"DP=(\d+);")
INDEL_RE = re.compile(r"INDEL")

def set_defaults():
  if len(sys.argv) == 2 and re.search(r'wgs|wes', sys.argv[1]):
    if sys.argv[1] == "wes":
      MAX_COV  = 200
      MIN_COV  = 10
    else: # wgs
      MAX_COV  = 100
      MIN_COV  = 5
  else:
    sys.stderr.write("Usage: " + sys.argv[0] + " wgs or wes")
    sys.exit(1)

def passes_filter(line):
  sl = line.split()

  # Qual filter
  col_qual = 5
  if float(sl[col_qual]) < MIN_QUAL:
    return False

  # Coverage
  match = COV_RE.search(line)
  if match:
    dp = int(match.group(1))
    if dp < MIN_COV and dp > MAX_COV:
      return False

  # Min number of samples having the snp
  match = SF_RE.search(line)
  if match:
    num_sfs = len(match.group(1).split(","))
    if num_sfs == 1:
      return False

  return True

def valid_snps():
  for line in sys.stdin:
    if line[0] == '#' or INDEL_RE.search(line): continue
    if passes_filter(line):
      yield(line)

def func_cons(line):
  match = EFF_RE.search(line)
  if match:
    return match.group(1)

def main():
  for line in valid_snps():
    print line.rstrip('\n')

if __name__=='__main__':
  set_defaults()
  sys.exit(main())
