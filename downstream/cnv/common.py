MIN_MAPQ = 10

def check_sorted(current_id, current_coor, prev):
  same_chrm = current_id == prev[0]
  prev_coor_bigger_than_current = current_coor < prev[1]
  if same_chrm and prev_coor_bigger_than_current:
    raise Exception('The bam does not seem to be coordinate sorted.')
  return (current_id, current_coor)

def good_alignment(ar):
  return not ar.is_unmapped and ar.mapq >= MIN_MAPQ and not ar.is_duplicate

def gen_chrm_lenghts(samfile):
  """ Given a bam (samfile), give me list of chrms and the sizes of chrms"""
  d_genome = {}
  chrm_names = [e['SN'] for e in samfile.header['SQ']]
  for name, size in zip(chrm_names, samfile.lengths):
    d_genome[name] = size
  return chrm_names, d_genome

def reads_at_locus(pucol, reads):
 for pileupread in pucol.pileups:
   al = pileupread.alignment
   if al.mapq >= 10 and not al.is_duplicate:
    reads.add(al.qname)
 return reads
