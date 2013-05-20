#!/usr/bin/env python
#
from __future__ import division
import sys
#
import drdcommon
from drdvcf import Vcf, VcfSnp
import drdplots as drdplots
from drdmath import std_dev
from groups import Group

usage = """
tool = snps_only_in_prog_or_cont.py

Given a merged vcf file and groups of samples, report snps where
the standard deviation of the variant allele ratios for all the groups
is over the min_threshold.

Usage:
  $ cat vcf | tool <groups_file> <min_std_threshold> <min_num_alleles_per_group> > out.vcf 2> stats.csv

groups.csv format:
        sample_id, group_name
        sp1, gr1
        sp2, gr2
        sp3, gr1
        ...
"""

class ComputeSnps():
  """Compute each of the snps
     STDIN: snp records that pass the filter.
     STDERR: var allele ratio per group and the final std_dev.
  """
  def __init__(self, vcf, o_groups, min_threshold, min_num_all):
    self.vcf, self.o_groups              = vcf, o_groups
    # min value for std and min number of calls per group to consider the snp
    self.min_threshold, self.min_num_all = min_threshold, min_num_all

  def run(self):
    self.dump_start()
    self.process_snps()

  def dump_start(self):
    print self.vcf.get_meta()      # dump the metadata + header from vcf
    for g in self.o_groups.groups: # csv output header
      sys.stderr.write("var_allele_ratio_group_%s " % g)
    sys.stderr.write("std_dev\n")

  def process_snps(self):
    """Per each snp, we want to decide if it is interesting. By interesting
    we mean that the samples in any group have a different genotype compare
    to the other groups.
    To do that we count the var allele freq per each of the samples in
    the group (#g1, #g2 ...) and apply the following condition:
    std_dev(#g1, #g2, #g3 ...) > X
    """
    ict = self.vcf.col_to_id
    for l in self.vcf.each_snp():
      self.h = self.o_groups.fresh_hash() # key: group , value: num of alt alleles seen
      snp = VcfSnp(l)
      self.process_genotypes(snp)
      self.check_filters_and_report(l)

  def process_genotypes(self, snp):
    for column, gts in snp.gtypes().items(): # column number and the tuple of genotypes per each sample
      for g in gts:                          # alleles in the current genotype
        sample_id    = self.vcf.col_to_id[column]
        sample_group = self.o_groups.what_is(sample_id)
        self.h[sample_group]["total"] += 1
        if g == 1: # We are seeing the var allele
          if sample_group:
            self.h[sample_group]["var_all"] += 1
          else:
            raise Exception('The sample (%s) does not belong to any group!' % sample_id)

  def check_filters_and_report(self, l):
    """Let's see if the snp passes the filters. If so, report it.
    """
    self.calculate_ratios(l)
    if self.enough_alleles(l) and self.std >= self.min_threshold:
      sys.stderr.write("%.3f %.3f %.3f\n" % (self.ratios[0], self.ratios[1], self.std))
      print l
      #print "\n%s\n%s\n" % (l, self.h)

  def enough_alleles(self, l):
    """We have to have a min number of alleles seen per each group (enough calls
    per each group)
    """
    for gn, n in self.h.items():
      if n["total"] < self.min_num_all: return False
    return True

  def calculate_ratios(self, l):
    """Compute the var allele ratios per each group and the std amoung all those ratios.
    """
    self.ratios = [] # ratio var / total
    for gn, n in self.h.items():
      self.ratios.append(n["var_all"]/n["total"])
    self.std = std_dev(self.ratios)

def do_work(fd_vcf, csv):
  vcf                = Vcf(fd_vcf)
  vcf.load_meta_header()
  grps               = Group(csv)
  #sys.stderr.write("# Of groups loaded: %s\n" % grps.num())
  min_threshold      = float(sys.argv[2])
  min_num_all        = float(sys.argv[3])
  ComputeSnps(vcf, grps, min_threshold, min_num_all).run()

def main():
  if len(sys.argv) == 4:
    fd_vcf = drdcommon.xopen("-")
    fd_csv = drdcommon.xopen(sys.argv[1])
    do_work(fd_vcf, fd_csv)
    fd_vcf.close()
    fd_csv.close()
  else:
    drdcommon.error("Incorrect # of params.", usage)

if __name__ == "__main__":
  main()
