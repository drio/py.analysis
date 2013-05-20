#!/usr/bin/env python
#
import json, requests
from drdcommon import *
from drdvcf import Vcf, VcfSnp

SERVER_URL = 'http://localhost:8080'

def doHeader(fd):
  vcf = Vcf(fd)
  vcf.load_meta_header()
  ot_already_there = str(vcf.check_info('OT'))
  if ot_already_there == True:
    error("This vcf seems to have an OT INFO already. Bailing out.")
  vcf.add_info('OT', '0', 'Flag', 'The site is on target.')
  # print "Checking if INFO id=OT is there: " + str(vcf.check_info('OT'))
  print vcf.get_meta()
  return vcf

def check_targets(chunk):
  # Check if the snps are on target or not
  data = { 'sites': [] }
  for l in chunk:
    chrm, coor = l.split()[0:2]
    data['sites'].append({'Chrm': chrm, 'Start': int(coor)})
  data_json = json.dumps(data)
  r = requests.post(SERVER_URL, data=data_json)

  # Add the OT INFO field if they are on target
  for i, ont in enumerate(json.loads(r.text)):
    if ont == 0: # not on target
      print chunk[i]
    else:
      print VcfSnp(chunk[i]).add_info('OT')

def doData(fd, vcf):
  chunk = []
  i = 0
  for l in vcf.each_snp():
    i+=1
    chunk.append(l)
    if i % 100000 == 0:
      check_targets(chunk)
      chunk = []
  check_targets(chunk) # We may have some left ...

def main():
  if len(sys.argv) == 2:
    fd = xopen(sys.argv[1])
    vcf = doHeader(fd)
    doData(fd, vcf)
    fd.close()
  else:
    error("Need vcf file.")

if __name__ == "__main__":
  main()
