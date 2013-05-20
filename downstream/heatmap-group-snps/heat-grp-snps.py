#!/usr/bin/env python
#
from __future__ import division
import pandas as pd
import numpy as np

from collections import defaultdict
import sys, os, re
#
import drdcommon
from drdvcf import Vcf, VcfSnp
import drdplots as drdplots
from drdmath import std_dev
from groups import Group

tool_name = os.path.basename(__file__)
usage = """
tool = %s

Given a VCF file, create a heatmap where each column is a SNP and
each row is the genotype for a particular sample. The posible
genotypes are: HOMO, HETE, NO_CALL, NO_COVERAGE

Usage:
  $ cat vcf | %s groups.tsv haplotype.tsv

groups.csv format:
        sample_id, group_name
        sp1, gr1
        sp2, gr2
        sp3, gr1
        ...
""" % (tool_name, tool_name)

def make_the_call(gts, _index, snp):
  """Given the dict of genotypes (gts) and the column index of the
  sample (_index) make a call:
  """
  we_have_a_call = gts.has_key(_index)
  HETE, HOMO, OTHER, NO_COV, NO_CALL = 0, 1, 2, 3, 4
  if we_have_a_call:
    a1, a2    = gts[_index][0], gts[_index][1]
    is_a_homo = (a1==1 or a2==1) and (a1==0 or a2==0)
    is_a_hete = a1==1 and a2==1
    if is_a_hete:
      return HETE
    elif is_a_homo:
      return HOMO
    else:
      return OTHER
  else:
    if snp.is_there_enough_coverage(_index):
      return NO_CALL
    else:
      return NO_COV

def prepare(vcf, grps_pheno, grps_haplo):
  """Prepare the data in the snps for the heatmap"""
  matrix, a_sites, a_groups = [], [], []

  for curr_grp in grps_pheno.groups:
    for _id in grps_pheno.indices_for_grp(curr_grp):
      pheno = curr_grp[0:5]
      haplo = grps_haplo.what_is(_id)
      a_groups.append(_id + "_" + pheno + "_" + haplo)

  for l in vcf.each_snp():
    snp     = VcfSnp(l)
    a_calls = []
    a_sites.append(snp.coordinate())
    gts     = snp.gtypes() # col_num -> gt_set
    for curr_grp in grps_pheno.groups:
      for _id in grps_pheno.indices_for_grp(curr_grp):
        _index = vcf.id_to_col[_id]
        a_calls.append(make_the_call(gts, _index, snp))
    matrix.append(a_calls)

  return np.transpose(np.array(matrix)), a_sites, a_groups

def do_work(fd_vcf, tsv_pheno, tsv_haplo):
  vcf  = Vcf(fd_vcf)
  vcf.load_meta_header()
  grps_pheno, grps_haplo = Group(tsv_pheno), Group(tsv_haplo)
  matrix, a_sites, a_groups = prepare(vcf, grps_pheno, grps_haplo)

  print matrix.shape
  print a_groups
  #values   = np.random.randn(100,100) * 10
  #a_sites  = ['1_100', '2_200', '3_300']
  #a_groups = ['grp1', 'grp1', 'grp2']
  cb_labels = ['HETE', 'HOMO_VAR', 'OTHER', 'NO_COVERAGE', 'HOMO_REF']
  drdplots.Heatmap(matrix, cb_labels, a_sites, a_groups).plot()

def main():
  if len(sys.argv) == 3:
    fd_vcf       = drdcommon.xopen("-")
    fd_pheno_tsv = drdcommon.xopen(sys.argv[1])
    fd_haplo_tsv = drdcommon.xopen(sys.argv[2])
    do_work(fd_vcf, fd_pheno_tsv, fd_haplo_tsv)
    fd_vcf.close()
    fd_pheno_tsv.close()
    fd_haplo_tsv.close()
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()

