#!/usr/bin/env python

import glob
import re
import sys
from drdcommon import xopen


class Data(object):
    def __init__(self):
        self.d = {}
        self.sep = "\t"

    def add(self, l, _id):
        d = self.d
        s = l.rstrip().split()
        #chrm    start   end     gene    exon_number     transcript_number       count   mean    std     min     25%     50%     75%     max
        chrm, start, end, gene, exon, trans, _, mean = s[0:8]
        chrm = re.sub("Chr", "", chrm)
        k = "%s %s %s %s %s %s" % (chrm, start, end, gene, exon, trans)
        if k not in d:
            d[k] = {}
        d[k][_id] = int(float(mean))

    def dump(self, ids):
        o = ""
        o += "chrm start end gene exon_number transcript_number "
        for i in ids:
            o += ("%s " % i)
        print re.sub("\s", self.sep, o)
        for k, means in self.d.items():
            o = k + " "
            for _id in ids:
                o += "%s " % means[_id]
            print re.sub("\s", self.sep, o)

d = Data()
ids = []
for f in glob.glob("*.stats.gz"):
    match = re.search("^(\d+)\.", f)
    if match:
        _id = int(match.group(1))
        first_line = True
        for l in xopen(f):
            if first_line:
                sys.stderr.write("%s\n" % _id)
                ids.append(_id)
                first_line = False
                continue
            else:
                d.add(l, _id)

d.dump(ids)
