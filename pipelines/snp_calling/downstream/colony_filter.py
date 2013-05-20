#!/usr/bin/env python

import csv, sys, re

def find_col_nums(row):
  cn, cn_id, cn_colony = 0, -1, -1
  for c in row:
    if c.upper() == 'DNAID': cn_id = cn
    if c.upper() == 'COLONY': cn_colony = cn
    cn+=1
  if cn_id == -1 or cn_colony == -1:
    err = 'I couldnt find the column numbers using the header. %d %d' % (cn_id, cn_colony)
    sys.stderr.write(err)
    raise
  return cn_id, cn_colony

def process_row(_id, colony, h):
  if colony not in h:
    h[colony] = set([])
  h[colony].add(_id)

def load_csv(fname):
  fd = open(fname, "r")
  reader = csv.reader(fd)
  i = 0
  h_colonies = {}
  for row in reader:
    if i == 0:
      cn_id, cn_colony = find_col_nums(row)
    else:
      process_row(row[cn_id], row[cn_colony], h_colonies)
    i+=1
  fd.close()
  return h_colonies

class SnpEntry:
  RE_ID = re.compile(r"(\d\/\d)")

  def __init__(self, line, itc):
    self.line = line
    self.id_to_gt = {}  # id to genotype
    self.itc      = itc # col_to_id from metadata in vcf
    self.process()

  def process(self):
    s = self.line.rstrip('\n').split('\t')
    self.chrm, self.pos = s[0], s[1]
	  # 0/1:10:139,0,7
    for i, gt in enumerate(s[9:]):
      match = self.RE_ID.search(gt)
      if match:
        self.id_to_gt[self.itc[i]] = match.group(1)
      else:
        self.id_to_gt[self.itc[i]] = '.'

  # Given a set of ids (colony) is the current snp present with the
  # same genotype in all the samples of the colony?
  def is_snp_in_colony_stringent(self, s_ids): # set of ids in colony
    n, gt = 0, "" # n_of_ids processed and genotype in colony
    for idx, _id in enumerate(s_ids):
      c_gt = self.id_to_gt[_id] # current genotype
      if c_gt == '.':
        return False
      if idx != 0 and c_gt != gt:
        return False
      gt = c_gt
      n  += 1
    if n == len(s_ids):
      return True
    else:
      return False

  # If the SNP is present in at least one animal, we consider that
  # a colony SNP
  # TODO: We are not looking at the actual genotype of the sample, we should
  # in order to decide if the snp belongs to the colony or not.
  #
  def is_snp_in_colony(self, s_ids): # set of ids in colony
    n, gt = 0, "" # n_of_ids processed and genotype in colony
    for idx, _id in enumerate(s_ids):
      c_gt = self.id_to_gt[_id] # current genotype
      if c_gt == '.':
        continue
      else:
        return True

class VCFile:
  RE_ID = re.compile(r"([\w_\d]+).bam$")

  def __init__(self, fd):
    self.fd        = fd # file pointer to input data (vcf)
    self.col_to_id = {} # what's the column number given a sample id ?
    self.load_header()

  def filter_by_colony(self, s_ids):
    for l in self.fd:
      snp = SnpEntry(l, self.col_to_id)
      if snp.is_snp_in_colony(s_ids):
        print l.rstrip('\n')
      #else: print snp.chrm, snp.id_to_gt, snp.is_snp_in_colony(s_ids)

  def load_header(self):
    for line in self.fd:
      line = line.rstrip('\n')
      if line[0] == '#':
        split = line.split()
        if split[0] == '#CHROM':
          sps = 9 # sample paths start column number
          self.sample_to_col_num(split[sps:]) # store the col number of the samples
          break # Next line will be a snp

  def sample_to_col_num(self, paths):
    for i, p in enumerate(paths): # index, path
      match = self.RE_ID.search(p)
      if match:
        _id = match.group(1)
        self.col_to_id[i] = _id
      else:
        raise 'Couldnt match id in header line.'

def usage(msg=""):
  if msg != "": sys.stderr.write("Ups! " + msg + "\n")
  sys.stderr.write("Usage: cat my.vcf | colony.py <id_to_colony.csv> <colony_name>")
  sys.stderr.write("\n")
  sys.exit(1)

def main():
  if len(sys.argv) != 3: usage('Wrong number of arguments: ' + str(sys.argv))
  csv_fname, colony_name = sys.argv[1], sys.argv[2]

  h_colonies = load_csv(csv_fname)
  if colony_name not in h_colonies:
    usage('I cannot find the colony: ' + colony_name)
  colony = h_colonies[colony_name]
  sys.stderr.write(colony_name + ": " + str(colony) + "\n")

  VCFile(sys.stdin).filter_by_colony(colony)

# Main
#
if __name__=='__main__':
  sys.exit(main())
