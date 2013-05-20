#!/usr/bin/env python
#
import json, requests
from drdcommon import *
from drdvcf import Vcf, VcfSnp

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

def load_cov(fd):
  cov = {}
  k = ""
  for l in fd:
    s = l.split()
    k = s[0]
    cov[k] = []
    for c in s[1:]: cov[k].append(c)
  return cov, len(cov[k])

def add_rdp(sl, a_cov):
  rdp = ";RDP=" + ",".join(a_cov)
  return "\t".join(sl[0:6] + [ sl[7].rstrip("\t") + rdp ] + sl[8:])

def do_data(vcf, h_cov, n_samples):
  for l in vcf.each_snp():
    s = l.split()
    k = s[0] + "_" + s[1]
    if k in h_cov:
      print add_rdp(s, h_cov[k])
    else:
      print add_rdp(s, ["0"] * n_samples)

def main():
  if len(sys.argv) == 2:
    sys.stderr.write("Loading coverage ...")
    h_cov, n_samples = load_cov(xopen('-'))
    sys.stderr.write("%d locations loaded.\n" % len(h_cov))
    fn    = sys.argv[1]
    fd    = xopen(fn)
    vcf   = doHeader(fd)
    do_data(vcf, h_cov, n_samples)
    fd.close()
  else:
    error("usage: cat cov.bed | add_rdp.py file.vcf")

if __name__ == "__main__":
  main()
