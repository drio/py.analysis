#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys
import os.path
import pysam
#
import drdcommon
from drdvcf import Vcf, VcfSnp

tool = os.path.basename(__file__)
usage = """
tool = %s

Given a vcf (stdin; snps only), compute the snp density

Usage:
  $ gzip -cd input.vcf.gz | grep -v "#" | %s > density.bed

""" % (tool, tool)


def do_work(fd_vcf):
    w_size = 10000
    w      = (0, 0+w_size)
    cc, pc = None, None
    num    = 0
    for l in fd_vcf:
        if l[0] != '#':
            cc, coor = l.split("\t")[0:2]
            coor = int(coor)
            if not pc:
                pc = cc
            if coor > w[1] or (pc and pc != cc):
                print "%s\t%s\t%s\t%s" % (pc, w[0], w[1], num)
                if pc != cc:
                    w = (0, 0+w_size)
                    if coor < w[1]:
                        num = 1
                    else:
                        num = 0
                else:
                    num = 1
                    w = (w[1]+1, w[1]+w_size)
            else:
                num += 1
            pc = cc
    print "%s\t%s\t%s\t%s" % (pc, w[0], w[1], num)


def main():
    if len(sys.argv) != 1:
        drdcommon.error("Wrong # of args", usage)
    if not drdcommon.data_in_stdin():
        drdcommon.error("No data in stdin.", usage)
    fd_vcf = drdcommon.xopen("-")
    do_work(fd_vcf)
    fd_vcf.close()

if __name__ == "__main__":
    main()
