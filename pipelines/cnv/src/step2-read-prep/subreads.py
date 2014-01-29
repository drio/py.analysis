#!/usr/bin/env python
#
# Given a sam stream, generate two subreads per each read
#

import sys
import re
import drdcommon

i = 1
for l in drdcommon.xopen("-"):
    seq = re.split(r'\t', l)[9]
    qua = re.split(r'\t', l)[10]
    _one_s, _two_s = seq[10:45+1], seq[46:81+1]
    _one_q, _two_q = qua[10:45+1], qua[46:81+1]
    sys.stdout.write("@%s\n%s\n+\n%s\n" % (i, _one_s, _one_q))
    sys.stdout.write(">%s\n%s\n+\n%s\n" % (i+1, _two_s, _two_q))
    i += 2
sys.stderr.write("%s sub-reads generated.\n" % i)
