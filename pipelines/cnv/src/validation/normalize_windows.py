#!/usr/bin/env python
#
# vim: ts=4 expandtab:
#
# normalize the windows regions from two data sets
#
# grep "Chr1" output_canavar.bed  | normalize_windows.py truth/bak/Dai.bed > normalized.chr1.bed
#
# Real example:
# $ cvar="filter/output.copynumber.bed"
# $ truth="truth/bak/*Dai*.bed"
# $ norm="/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/validation/normalize_windows.py"
# $ grep -P "^18" $cvar |  cut -f1,2,3,5 | awk '{print "chr"$_}' | $norm $truth | tail
#
import sys
import numpy as np
import drdcommon as drd

three_hundre_mil = 300000000
a = np.zeros(shape=(three_hundre_mil))

# Save the metric values for all the bp from the stdin input
working_chrm = ""
for i in drd.xopen("-"):
    working_chrm, start, stop, val = i.split()
    for j in range(int(start), int(stop)+1):
        a[j] = val

# Iterate over second output and generate new windows for first input
for i in drd.xopen(sys.argv[1]):
    chrm, start, stop, val = i.split()
    if working_chrm == chrm:
        s, e = int(start), int(stop)
        print "%s\t%s\t%s\t%s" % (chrm, start, stop, int(np.median(a[s:e+1])))

# (cd /Users/drio/dev/py.analysis; make update)

