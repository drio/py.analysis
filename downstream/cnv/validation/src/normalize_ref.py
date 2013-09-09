#!/usr/bin/env python
#
import sys, os.path
import drdcommon

max_c = int(sys.argv[1]) if len(sys.argv) == 2 else 80
nc = -1
first = True
for l in sys.stdin:
  if l[0] == '>':
    if first:
      first=False
    else:
      print ""
    sys.stdout.write(l)
    nc = 0
  else:
    for c in l.rstrip():
      nc += 1
      sys.stdout.write(c)
      if nc == max_c:
        print ""
        nc = 0
print ""
