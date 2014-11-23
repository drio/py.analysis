#!/usr/bin/env python
"""
Find gaps coordinates.
Example:
    $ pad_gaps fasta.fa > gaps_coor.bed
"""
import sys
import drdcommon as drd

chrm, coor, start = None, None, None
for l in drd.xopen(sys.argv[1]):

    l = l.strip()
    if l[0] == '>':
        if start:
            print "%s\t%s\t%s" % (chrm, start-1, coor-1)
            start = None
        chrm = l[1:]
        coor = 1
    else:
        for c in l:
            if not start and c == 'N':
                start = coor

            if start and c != 'N':
                print "%s\t%s\t%s" % (chrm, start, coor-1)
                start = None

            coor += 1

