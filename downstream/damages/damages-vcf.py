#!/usr/bin/env python
#
import sys
import drdcommon
import os.path

tool = os.path.basename(__file__)
usage = """
tool = %s

Usage:
  $ %s <sample_id> <vcf>

""" % (tool, tool)

def main():
    if len(sys.argv) == 3:
        sample_id, vcf_file = sys.argv[1:]

        d = {}
        for l in drdcommon.xopen(vcf_file):
            if l[0] == "#":
                continue
            v_chrm, v_coor = l.split("\t")[0:2]
            d[v_chrm + "_" + v_coor] = True

        for l in drdcommon.xopen("-"):
            l = l.rstrip()
            s = l.split("\t")
            chrm, coor = s[0:2]
            if (chrm + "_" + coor) in d:
                print sample_id + "\t" + l
    else:
        print usage

if __name__ == '__main__':
    main()
