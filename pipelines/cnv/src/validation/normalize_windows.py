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
import math
import numpy as np
import scipy
import drdcommon as drd

_250m = 250000000

def log(msg):
    sys.stderr.write(">> " + msg + "\n")

def method1(input1, input2):
    """
           4      3              5            1
    i1  |------|-------|------------------|--------|
              7                10               1
    i2  |----------|-------------------------|-----|
             3.5             4                  1
    out |----------|-------------------------|-----|
    """
    a = np.zeros(shape=(three_hundre_mil))

    # Save the metric values for all the bp from the stdin input
    working_chrm = ""
    for i in input1:
        working_chrm, start, stop, val = i.split()
        for j in range(int(start), int(stop)+1):
            a[j] = val

    # Iterate over second output and generate new windows for first input
    for i in drd.xopen(input2):
        chrm, start, stop, val = i.split()
        if working_chrm == chrm:
            s, e = int(start), int(stop)
            print "%s\t%s\t%s\t%s" % (chrm, start, stop, int(np.median(a[s:e+1])))

def compute_vals_wins(stream, working_chrm):
    """
    Represent in an array the windows' coordinates
    Save the metric values per each bp

    Input:
        chr1 0 4 5
        chr1 4 8 7
    Returns:
        012345678
        555557777 (a_wins)
    and:
        100010001 (a_vals)
    """
    a_vals = scipy.zeros((_250m), float)
    a_wins = scipy.zeros((_250m), int)
    for i in stream:
        chrm, start, stop, val = i.rstrip().split()
        if working_chrm == chrm:
            s, e = int(start), int(stop)
            a_wins[s], a_wins[e] = 1, 1
            for j in range(s, e+1):
                a_vals[j] = math.ceil(float(val)*100)/100
    return a_vals, a_wins

def cv(a, s, e):
    return a[s:e+1][1]

def method2(stream1, stream2, working_chrm):
    """
               5            10                 4
    s1   |-----------|--------------|---------------------|
            4                    2                  6
    s2   |------|-----------------------------|-----------|

    out  |------|----|--------------|---------|-----------|
           4,5    2,5      10,2         4,2       4,6
    """
    log("Creating arrays for first stream")
    a_s1_vals, a_s1_wins = compute_vals_wins(stream1, working_chrm)
    log("Creating arrays for second stream")
    a_s2_vals, a_s2_wins = compute_vals_wins(stream2, working_chrm)

    log("Finding coordinate locations")
    o_wins = a_s1_wins | a_s2_wins
    log("Iterating over %s windows" % len(o_wins))
    _first, _prev = True, None
    for coor in np.where(o_wins == 1)[0]:
        if _first:
            _prev = coor
            _first = False
        else:
            s, e = _prev, coor
            print "%s %s %s %s %s" % (working_chrm, s, e, cv(a_s1_vals, s, e), cv(a_s2_vals, s, e))
            _prev = coor

if __name__ == "__main__":
    # method1(drd.xopen(sys.argv[1]), drd.xopen(sys.argv[2]))
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: tool <bed_file1> <bed_file2> <chromosome>\n")
        exit(1)
    method2(drd.xopen(sys.argv[1]), drd.xopen(sys.argv[2]), sys.argv[3])

