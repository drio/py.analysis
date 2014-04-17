#!/usr/bin/env python
#
# vim: set ts=4 et:
#
# Use the counts of how many samples had high or
# low coverage to call the under/over representation
#
import sys
import drdcommon

def key(chrm, start, end):
    return "%s_%s_%s" % (chrm, start, end)

d = {}
for l in drdcommon.xopen("-"):
    chrm, start, stop, _min, _max, _pass = l.rstrip().split("\t")
    k = key(chrm, start, stop)

    d[k] = [int(_min), int(_max), int(_pass)]

for l in drdcommon.xopen(sys.argv[1]):
    l = l.rstrip()
    chrm, start, stop, _ = l.rstrip().split("\t")
    k = key(chrm, start, stop)
    if k in d:
        _min, _max, _pass = d[k]
        if _min > _max:
            print l + "\t" + "under"
        else:
            print l + "\t" + "over"
    else:
        sys.stderr.write("I cannot find key!")
        raise
