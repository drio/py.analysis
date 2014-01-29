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
    _one, _two = seq[10:45+1], seq[46:81+1]
    sys.stdout.write(">%s\n%s\n" % (i, _one))
    sys.stdout.write(">%s\n%s\n" % (i+1, _two))
    i += 2
sys.stderr.write("%s sub-reads generated, # of reads in input stream = %s\n" % (i, i/2))
