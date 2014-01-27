#!/bin/bash
#
# CNV pipeline: First step
# Masking over-represented parts of the genome
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

error() {
  local msg=$1
  echo $msg
  exit 1
}

fasta_genome=$1

# Some sanity checks
#
[ ".$fasta_genome" == "." ] && error "Need path to fasta file genome."
for b in R grep cut awk mrsfast maskFastaFromBed fastaFromBed
do
  `which $b &>/dev/null` || error "$b not in path."
done

# Gen chrm info
#################################################################################
_chrm_info_bed="chrm_info.bed"
if [ ! -f $_chrm_info_bed ]
then
  cat << _EOF
I can't find $_chrm_info_bed (chrm_name start end). You may find this useful:
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from rheMac2.chromInfo" | \
grep chr | awk '{print \$1"\t"\$2"\t"\$2}' | grep -v size > $_chrm_info_bed
_EOF
exit 1
fi

# Make index
#################################################################################
echo -e "mrsfast --index $fasta_genome --ws 12 ${fasta_genome}.index\tindexing\t-"
echo

# Generate kmer intervals (bed)
#################################################################################
K=36
STEP=5
while read line
do
	_cName=$(echo -e $line | cut -f1 -d' ')
	_cLen=$(echo -e $line | cut -f2 -d' ')
	echo -e "$SRC_DIR/makeIntervalsBED_v2.pl $_cName $_cLen $K $STEP > ${_cName}_k${K}_step${STEP}_intervals.bed\tintervals\t-";
done < $_chrm_info_bed
echo

# Kmerify reference
#################################################################################
while read line
do
	_cName=$(echo -e $line | cut -f1 -d' ')
  _bed="${_cName}_k${K}_step${STEP}_intervals.bed"
  _out="${_cName}.fa"
  echo -e "fastaFromBed -fi $fasta_genome -bed $_bed -fo $_out\tkmerify\tintervals"
done < $_chrm_info_bed
echo

# Map kmers against genome
#################################################################################
while read line
do
	_cName=$(echo -e $line | cut -f1 -d' ')
  _reads="${_cName}.fa"
  _out="${_cName}_k${K}_step${STEP}_intervals.map"
  echo -e "mrsfast --search $fasta_genome --seq $_reads -o $_out -e 2\tmapping\tkmerify"
done < $_chrm_info_bed
echo

# Count how many hits we have per kmer
#################################################################################
rm -f $_counts_bed
while read line
do
	_cName=$(echo -e $line | cut -f1 -d' ')
  _kmers="${_cName}.fa"
  _maps="${_cName}_k${K}_step${STEP}_intervals.map"
  _counts_bed=counts._${_cName}.bed
  echo -e "${SRC_DIR}/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl $_kmers $_maps $_cName > $_counts_bed\tcounts\tintervals,mappings"
done < $_chrm_info_bed
echo

# Count how many hits we have per kmer
#################################################################################
_merged_counts_bed="counts.bed"
echo -e "cat counts* | grep -v \# | grep -v -P \^\$ >> $_merged_counts_bed\tmerge_counts\tcounts"
echo

#cat $fasta_genome | src/step1/calculateAssemblyStats_v2.pl > stats.txt

# Mask genome
#################################################################################
echo -e "maskFastaFromBed -fi $fasta_genome -bed $_merged_counts_bed -fo masked.fa\tmasking\tmerge_counts"


# Make dot plot (cumulative percentage of the genome covered for different numbers
# of hits allowed per each kmer.
#################################################################################
echo -e "cat $_merged_counts_bed | R CMD BATCH ${SRC_DIR}/makeHistogramCounts.R\tdotplot\tmerge_counts"
echo
