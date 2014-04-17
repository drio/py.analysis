#!/usr/bin/env python
#
# Given the enumeration from enumerate.sh as input,
# report, per each location, how many samples we have
# each filtering category
#
import sys
import drdcommon

d = {}
_a = {}
for l in drdcommon.xopen("-"):
  _id, chrm, start, end, _type = l.rstrip().split("\t")

  if _id not in _a:
    sys.stderr.write(_id + "\n")
    _a[_id] = {}

  k = "%s_%s_%s" % (chrm, start, end)

  if k not in _a[_id]:
    _a[_id][k] = True

    if not k in d:
      d[k] = {}
      for _t in [ "min", "max", "pass" ]:
        d[k][_t] = 0

    d[k][_type] += 1

print "chrm start stop min max pass".replace("\s", "\t")
for k, v in d.items():
  s = k.split("_")
  if len(s) != 3:
    sys.stderr.write("Error: k = %s\n" % k)
  else:
    chrm, start, end = s
    _ = "%s %s %s %s %s %s" % (chrm, start, end, v["min"], v["max"], v["pass"])
    print _.replace(" ", "\t")

