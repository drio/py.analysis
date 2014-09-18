#!/usr/bin/env python
#
import sys
import drdcommon
import os.path

tool = os.path.basename(__file__)
usage = """
tool = %s

Usage:
  $ %s <lift_over_file>

Format:

chrm start end chrm start end
""" % (tool, tool)

def main():
    if len(sys.argv) == 2:
        lo_file = sys.argv[1]
        link = {}
        for l in drdcommon.xopen(lo_file):
            s = l.rstrip().split("\t")
            chrm_from, start_from, end_from = s[0:3] # hsap
            chrm_to, start_to, end_to = s[3:6] # rhmac
            # human -> rhmac
            #link[chrm_from + "_" + end_from] = chrm_to + "_" + end_to
            link[chrm_from + "_" + start_from] = chrm_to + "_" + start_to

        for l in drdcommon.xopen("-"):
            s = l.rstrip().split("\t")
            key_hsap = "_".join(s[0:2])
            if key_hsap in link:
                rh_coor = link[key_hsap]
            else:
                rh_coor = "-_-"
            print rh_coor.replace("_", "\t") + "\t" + l.rstrip()
    else:
        print usage

if __name__ == '__main__':
    main()
