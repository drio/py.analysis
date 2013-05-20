#!/usr/bin/env python
#
import json, requests, logging, os
from drdcommon import *
from drdvcf import Vcf, VcfSnp

tool = os.path.basename(__file__)
logging.basicConfig(format=tool + ': %(asctime)s >> %(message)s', level=logging.DEBUG)

def doHeader(fd):
  vcf = Vcf(fd)
  vcf.load_meta_header()
  ot_already_there = str(vcf.check_info('RDP'))
  if ot_already_there == True:
    error("This vcf seems to have an RDP INFO already. Bailing out.")
  vcf.add_info('RDP', '1', 'Integer', 'Raw read coverage at locus.')
  # print "Checking if INFO id=OT is there: " + str(vcf.check_info('OT'))
  print vcf.get_meta()
  return vcf

def add_rdp(sl, a_cov):
  """ sl: the split line for the snp
      a_cov: is the array of depth of coverage
  """
  rdp = ";RDP=" + ",".join(a_cov)
  return "\t".join(sl[0:6] + [ sl[7].rstrip("\t") + rdp ] + sl[8:])

def pop_coor(fd_cov):
  s_cov       = fd_cov.readline().split()
  coor, a_cov = s_cov[0], s_cov[1:]
  return coor, a_cov

def raise_coor_do_not_match(coor_cov, snp_coor):
  msg = '''Coordinates in coverage file do not match coor. in vcf file.
        coor_cov(%s) != snp(%s) ''' % (coor_cov, snp_coor)
  raise Exception(msg)

def print_snps(a_snps, coor_cov, a_cov, snp_line):
  for s in a_snps:
    if coor_cov != s.coordinate():
      raise_coor_do_not_match(coor_cov, s.coordinate())

def process_snps(vcf, fd_cov):
  """ We have to read each coverage line and the corresponding snp.
  """
  for l in vcf.each_snp():
    snp = VcfSnp(l)
    coor_cov, a_cov = pop_coor(fd_cov)
    print add_rdp(l.split(), a_cov)

def main():
  if len(sys.argv) == 2 and data_in_stdin():
    fd_vcf  = xopen(sys.argv[1])
    fd_cov  = xopen("-")
    vcf     = doHeader(fd_vcf)
    process_snps(vcf, fd_cov)

    fd_cov.close()
    fd_vcf.close()
    logging.info("Done.")
  else:
    error("usage: cat cov.bed | add_rdp.py file.vcf")

if __name__ == "__main__":
  main()
