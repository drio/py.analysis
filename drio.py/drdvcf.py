# ##fileformat=VCFv4.1
# ##fileDate=20090805
# ##source=myImputationProgramV3.1
# ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
# ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
# ##phasing=partial
# ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
# ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
# ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
# ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
# ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
# ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
# ##FILTER=<ID=q10,Description="Quality below 10">
# ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
# ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
# #CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
# 20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
# 20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3

#Chr1    424     .       A       T       85.00   .       AC1=1;AC=1;AF1=0.5;AN=2;DP4=1,1,2,2;DP=6;EFF=INTERGENIC(MODIFIER|||||||||);FQ=31;MQ=60;PV4=1,0.32,1,1;SF=24;VDB=0.0322  GT:GQ:PL        .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .    ..       .       .       .       .       .       .       0/1:61:115,0,58 .       .       .       .       .       .       .       .       .       .       .       .       .       .       .

import re

VCF_SNP_AN  = re.compile(r"AN=(\d+)")       # match AN
VCF_SNP_AC  = re.compile(r"AC=([\d,]+)")    # match AC
VCF_SNP_DP4 = re.compile(r"DP4=([\d,]+)")   # match DP4
SAMPLE_GTS  = re.compile(r"^(\d+)/(\d+):")  # match for sample genotypes
SNP_EFF     = re.compile(r"EFF=")           # match for snp effect
FUNC_CONS   = re.compile(r"EFF=(\w+)\(")    # match for the functional consequence
RDP         = re.compile(r"RDP=([\d,]+)")   # match RDP
NON_CODING  = re.compile(r"UTR|NONE|INTERGENIC|INTRON|UPSTREAM|DOWNSTREAM") # non coding region
INDEL       = re.compile(r"INDEL")

class VcfSnp:
  INFO_COLUMN = 7
  GT_COLUMN   = 9 # genotypes column

  def __init__(self, line):
    self.line = line
    self.s_line = line.split()
    self.set_anno_info()
    self.load_rdp()

  def species_snp(self):
    """
    If we align against a more distant reference than the genome
    of the dataset, we will have a high number of snps that are
    specific of that species. We may want to remove those.
    For example, mapping a cyno against rhesus macaque.
    This method tells us if self is one of those type of snps. We just
    have to check if all the samples have the same snp/genotype.
    """
    gtypes = self.gtypes()
    try:
      first_gt_set = gtypes.items()[0][1]
    except:
      print self.line
      print gtypes
      raise
    for col, c_gt in gtypes.items():
      if c_gt != first_gt_set:
        return False
    return True

  def load_rdp(self):
    match = RDP.search(self.line)
    if match:
      hits = match.group(1)
      self.rdp = [int(x) for x in hits.split(",")]
    else:
      self.rdp = []

  def is_an_indel(self):
    match = INDEL.search(self.line)
    return True if match else False

  def is_a_substitution(self):
    return not self.is_an_indel()

  def has_low_quality(self, max_qual):
    return float(self.s_line[5]) < max_qual

  def has_high_quality(self, min_qual):
    return float(self.s_line[5]) > min_qual

  def in_coding_region(self):
    match = NON_CODING.search(self.line)
    return False if match else True

  def coordinate(self, delim="_"):
    return self.s_line[0] + delim + self.s_line[1]

  def add_info(self, info_str):
    cols      = self.line.split()
    new_infos = cols[self.INFO_COLUMN] + ";" + info_str
    new_line  = cols[0:self.INFO_COLUMN] + [new_infos] + cols[self.INFO_COLUMN+1:-1]
    self.line = "\t".join(new_line)
    return self

  def alternative_allele_counts(self):
    match = VCF_SNP_AC.search(self.line)
    if match:
      counts = match.group(1)
      return [int(x) for x in counts.split(",")]

  def total_num_alleles(self):
    match = VCF_SNP_AN.search(self.line)
    if match:
      an = match.group(1)
      return int(an)

  def ref_var_num_alleles(self):
    match = VCF_SNP_DP4.search(self.line)
    if match:
      counts = match.group(1)
      refp, refn, varp, varn = [int(x) for x in counts.split(",")]
      return refp, refn, varp, varn

  def dp4_var_allele_ratio(self):
    match = VCF_SNP_DP4.search(self.line)
    if match:
      counts = match.group(1)
      refp, refn, varp, varn = [int(x) for x in counts.split(",")]
      return round((0. + varp + varn) / (refp + refn + varp + varn), 2)
    else:
      raise('DP4 not found in SNP data.')

  def is_there_enough_coverage(self, col_num, min_num_reads=4):
    if self.we_have_rdp():
      return self.rdp[col_num] >= min_num_reads
    else:
      #from random import randint
      #randint(0,100)
      raise Exception('RDP not available for this dataset.')

  def we_have_rdp(self):
    return len(self.rdp) > 0

  def gtypes(self):
    """return an array with the genotypes per each sample
       gts[column_number] = (0,1)
       column_number of the sample
    """
    gts = {}
    for col_num, cv in enumerate(self.line.split()[self.GT_COLUMN:]):
      match = SAMPLE_GTS.search(cv)
      if match:
        a1, a2 = match.group(1), match.group(2)
        gts[col_num] = (int(a1), int(a2))
      else:
        if cv != '.':
          raise Exception('VcfSnp.gtypes(): Unexpected entry in sample genotype: [%s]' % cv)
    return gts

  def set_anno_info(self):
    self.annotated, self.func_cons, self.gene, self.impact = [ False, "", "", "" ]
    for c in self.line.split(";"):
      match = SNP_EFF.search(c)
      if match: # c contains the whole INFO column
        self.annotated = True
        self.effect = c

        # func_cons
        match = FUNC_CONS.search(c)
        if match: self.func_cons = match.group(1)

        # gene
        col_gene_location = 5
        split_col = c.split("|")
        self.gene = split_col[col_gene_location]

        # impact
        col_impact_location = 0
        self.impact = split_col[col_impact_location].split("(")[1]

  def __str__(self):
    return self.line

class Vcf:
  RE_INFO_ID = re.compile(r"##INFO=<ID=(\w+),") # match info id
  INFO_SEED  = '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">'

  def __init__(self, fd): # file descriptor to the vcf file
    self.fd     = fd
    self.meta   = []
    self.header = []

  def get_meta(self):
    return '\n'.join(self.meta + [self.header])

  def load_meta_header(self):
    for line in self.fd:
      line = line.rstrip('\n')
      if line[0:2] == '##': #metadat
        self.meta.append(line)
      else: # header
        self.header = line
        # Save the col number per each sample
        self.col_to_id = {}
        self.id_to_col = {}
        self.num_of_samples = 0
        for index, _id in enumerate(line.split()[9:]):
          self.col_to_id[index] = _id
          self.id_to_col[_id]   = index
          self.num_of_samples += 1
        break

  def check_info(self, _id): # Do we have an INFO entry with this id ?
    for l in self.meta:
      match = self.RE_INFO_ID.search(l)
      if match:
        m_id = match.group(1)
        if m_id == _id: return True
    return False

  def add_info(self, _id, num, _type, desc): # Add a new INFO to metadata
    if self.check_info(_id):
      raise "INFO with id: %s is already in the metadata" % _id
    else:
      new_info = self.INFO_SEED % (_id, num, _type, desc)
      self.meta.append(new_info)

  def each_snp(self):
    for line in self.fd:
      if line[0] != '#':
        yield line.rstrip('\n')
















