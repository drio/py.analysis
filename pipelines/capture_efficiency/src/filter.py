#!/usr/bin/env python

import sys
from drdcommon import xopen

first_line = True
d = []
for l in xopen("-"):
    if first_line:
        first_line = False
        continue
    else:
        d.append((float(l.split()[7]), l.rstrip()))

_sum = 0
for t in d:
    _sum += t[0]
mean = _sum/len(d)

_min, _max = mean-(mean-20), mean+(mean-20)
drops = 0
for t in d:
    m, line = t
    if m < _min or m > _max:
        drops += 1
        _ = ""
        if m < _min:
            _ = "min"
        else:
            _ = "max"
        sys.stderr.write(line + "\t" + _ + "\n")
    else:
        print line + "\t" + "pass"
