#!/usr/bin/env python
#
from __future__ import division
from collections import defaultdict
import sys, os, re
import os.path
import logging
from optparse import OptionParser

import drdcommon
from drdvcf import Vcf, VcfSnp
from drdmath import std_dev

tool_name = os.path.basename(__file__)
usage = """
Given a vcf stream (stdin) compute the snp frequency.

  $ %s <-wes|-wgs> -v <vcf_file> [-s] [-l]

Notes:
  For wgs we will use genome size of: 3Gbp
  and 34Mbp for wes.""" % (tool_name)

def raise_it(msg):
  raise Exception(msg)

class SnpFreq(object):
  GENOME_SIZE = {'wgs': 3000000000, 'wes': 34000000}
  MIN_QUAL    = 20
  def __init__(self, fd_vcf, exp_type, options):
    self.fd_vcf   = fd_vcf
    self.exp_type = exp_type
    self.drop = options.drop
    self.list_s_snps = options.list_s_snps
    self.__validate_type()
    self.coordinates_in_file = options.coordinates_in_file
    self.__load_vcf()
    if options.coor_fn:
      self.coor_fn = options.coor_fn
      self.__load_species_snp_coordinates()

  def __load_species_snp_coordinates(self):
    fd = drdcommon.xopen(self.coor_fn)
    d = {}
    self.d_species_coor = d
    n = 0
    for l in fd:
      n += 1
      chrm, coor = l.split()
      if not d.has_key(chrm):
        d[chrm] = {}
      d[chrm][int(coor)] = 1
    fd.close()
    logging.info("# of coordinates loaded: %d" % n)
    logging.info("current memory usage in %dkb" % drdcommon.memory_usage())

  def __load_vcf(self):
    self.vcf = Vcf(self.fd_vcf)
    self.vcf.load_meta_header()

    if self.drop and (not self.coordinates_in_file and self.vcf.num_of_samples < 2):
      drdcommon.error("I need a population level vcf in order to drop species snps.")

  def __validate_type(self):
    if self.exp_type != 'wgs' and \
       self.exp_type != 'wes' and \
       self.exp_type != 'null':
      raise_it('Invalid experiment type. Valid types: wgs or wes')

  def __list_species_snps(self):
    for l in self.vcf.each_snp():
      snp = VcfSnp(l)
      if snp.is_a_substitution() and \
         snp.has_high_quality(self.MIN_QUAL) and \
         snp.species_snp():
        print(snp.coordinate(' '))

  def __is_a_species_snp(self, snp):
    if self.coordinates_in_file:
      ch, co = snp.coordinate(' ').split()
      return self.d_species_coor.has_key(ch) and self.d_species_coor[ch].has_key(int(co))
    else:
      return snp.species_snp()

  def __calculate_snp_freq(self):
    """
    Compute the snp frequency (# of snps per kbp)
    Drop snps that are indels, have low quality
    If wes, also drop non coding regions
    If drop is True, we have to drop species snps
    """
    num_snps = 0
    total = 0
    for l in self.vcf.each_snp():
      snp = VcfSnp(l)
      total += 1
      if snp.is_a_substitution() and snp.has_high_quality(self.MIN_QUAL):
        if self.exp_type == 'wgs':
          if not self.drop or (self.drop and not self.__is_a_species_snp(snp)):
            num_snps += 1
        if self.exp_type == 'wes' and snp.in_coding_region():
          if not self.drop or (self.drop and not self.__is_a_species_snp(snp)):
            num_snps += 1

    logging.info("Total/counted: %d/%d" % (total, num_snps))
    return (float(num_snps)/self.GENOME_SIZE[self.exp_type])*1000

  def run(self):
    if self.list_s_snps:
      self.__list_species_snps()
    else:
      return self.__calculate_snp_freq()

class ProcessArgs(object):
  def __init__(self):
    self.process()
    self.check_logic()

  def process(self):
    parser = OptionParser(usage)

    parser.add_option("-g", "--wgs",
      action="store_true", dest="wgs", default=False,
      help="We are working with wgs data.")

    parser.add_option("-e", "--wes",
      action="store_true", dest="wes", default=False,
      help="We are working with wes data.")

    parser.add_option("-s", "--drop_species",
      action="store_true", dest="drop", default=False,
      help="Do not include species snps in the calculations.")

    parser.add_option("-v", "--vcf", dest="vcf_fn",
      help="vcf file_name. Use '-' for stdin.")

    parser.add_option("-l", "--list_s_snps",
      action="store_true", dest="list_s_snps", default=False,
      help="List species snps in stdin")

    parser.add_option("-c", "--coordinates_drop", dest="coor_fn",
      help="coordinates to drop file_name. Use '-' for stdin.")

    self.options, self.args = parser.parse_args()
    self.parser = parser

  def check_logic(self):
    options, args = self.options, self.args

    if not options.vcf_fn:
      drdcommon.error("Need vcf file.")

    if options.vcf_fn == '-' and not drdcommon.data_in_stdin():
      drdcommon.error("No data in stdin.")

    if not options.vcf_fn == '-' and not os.path.isfile(options.vcf_fn):
      drdcommon.error("Vcf file does not exists.")

    if options.coor_fn:
      if options.coor_fn == '-' and options.vcf_fn == '-':
        drdcommon.error("I cannot read two streams from stdin.")
      if not os.path.isfile(options.coor_fn):
        drdcommon.error("coor file does not exists.")
      options.drop = True
      options.coordinates_in_file = True
    else:
      options.coordinates_in_file = False

    if len(args) == 0:
      if options.wes:
        self.exp_type = 'wes'
      elif options.wgs:
        self.exp_type = 'wgs'
      else:
        self.exp_type = 'null'
        if not self.options.list_s_snps:
          drdcommon.error("Experiment type not set.")
    else:
      drdcommon.error("Incorrect # of params.")

  def run(self):
    logging.basicConfig(level=logging.INFO)
    fd_vcf = drdcommon.xopen(self.options.vcf_fn)
    sf = SnpFreq(fd_vcf, self.exp_type, self.options)
    if self.options.list_s_snps:
      sf.run()
    else:
      print "%f" % sf.run()
    fd_vcf.close()

def main():
  ProcessArgs().run()

if __name__ == "__main__":
  main()
