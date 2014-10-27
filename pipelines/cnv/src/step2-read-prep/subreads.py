#!/usr/bin/env python
#
# Given a fq/sam stream, generate two subreads per each read
#

import sys
import re
import drdcommon

def gen_subreads(seq, qua):
    #seq = re.split(r'\t', l)[9]
    #qua = re.split(r'\t', l)[10]
    return [ (seq[10:45+1], qua[10:45+1]), (seq[46:81+1],qua[46:81+1]) ]


def dump(i, sub_reads):
    for idx, t in enumerate(sub_reads):
        s, q = t
        sys.stdout.write("@%s\n%s\n+\n%s\n" % (idx+i, s, q))


i, sr = 1, 0
seq, qual = None, None
for l in drdcommon.xopen("-"):
    if sr == 1:
        seq = l.rstrip()

    if sr == 3:
        qual = l.rstrip()
        dump(i, gen_subreads(seq, qual))
        i += 2
        sr = -1

    sr += 1

sys.stderr.write("%s sub-reads generated.\n" % i)
