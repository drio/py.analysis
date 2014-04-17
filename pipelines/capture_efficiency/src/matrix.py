#!/usr/bin/env python

import glob
import re
import sys
from drdcommon import xopen


class Data(object):
    def __init__(self):
        self.d = {}
        self.sep = "\t"

    def add(self, l, _passes):
        d = self.d
        s = l.rstrip().split()
        #chrm    start   end     gene    exon_number     transcript_number       count   mean    std     min     25%     50%     75%     max
        chrm, start, end, gene, exon, trans, _, mean = s[0:8]
        chrm = re.sub("Chr", "", chrm)
        k = "%s %s %s %s %s %s" % (chrm, start, end, gene, exon, trans)
        if k not in d:
            d[k] = 0
        d[k] += _passes

    def dump(self):
        o = ""
        o += "chrm start end gene exon_number transcript_number num_samples_pass"
        print re.sub("\s", self.sep, o)
        for k, num_p in self.d.items():
            o = k + " " + str(num_p)
            print re.sub("\s", self.sep, o)


def process_file(wild_file, d, _passes):
    for f in glob.glob(wild_file):
        match = re.search("^(\d+)\.", f)
        if match:
            _id = int(match.group(1))
            first_line = True
            for l in xopen(f):
                if first_line:
                    sys.stderr.write("%s\n" % _id)
                    first_line = False
                    continue
                else:
                    d.add(l, _passes)
    return d

d = Data()
for f in [("*pass*.gz", 1), ("*out*.gz", 0)]:
    d = process_file(f[0], d, f[1])
d.dump()
