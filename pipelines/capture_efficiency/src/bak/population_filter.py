#!/usr/bin/env python

import sys
from drdcommon import xopen

min_num_samples = int(sys.argv[1])

first_line = True
for l in xopen("-"):
    if first_line:
        first_line = False
        continue
    else:
        #chrm    start   end     gene    exon_number     transcript_number       32510 ..
        num = 0
        for i in l.strip().split("\t")[6:]:
            i = int(i)
            if i > 0:
                num += 1


        out = "\t".join(l.strip().split("\t")[0:6]) + "\t" + str(num) + "\n"

        if num >= min_num_samples:
            sys.stdout.write(out)
        else:
            sys.stderr.write(out)
