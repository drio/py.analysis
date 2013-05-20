#!/usr/bin/env python
#
# show how many snps we have per each func cons type
#
import sys, re

EFF_RE   = re.compile(r"EFF=([\w_]+)\(")

def func_cons(line):
  match = EFF_RE.search(line)
  if match:
    return match.group(1)

def main():
  dict_func = {}
  total = 0
  for line in sys.stdin:
    if line[0] == '#': continue
    cons = func_cons(line)
    if cons == None:
      print "ERROR: I cannot extract func_cons for: "
      print line
      os.exit(1)
    if dict_func.has_key(cons):
      dict_func[cons]+=1
    else:
      dict_func[cons]=1
    total+=1

  dict_func["TOTAL"] = total
  for key, val in dict_func.iteritems():
    print key, val

if __name__=='__main__':
  sys.exit(main())
