#!/bin/bash
#
[ ! -f "../common.sh" ] && echo "I can't find ../common.sh" && exit 1

src_dir=$(dirname $(readlink -f $0))
source ../common.sh

mkdir -p steps

index() {
  local f_to_index=$1
  mrfast --index $f_to_index
}

gen_chrm_info() {
  local genome="hg19"
  mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from $genome.chromInfo" | \
  grep chr | awk '{print $1"\t"$2"\t"$2}' | grep -v size > $chrm_info
}

mask_rtf_simple_gaps() {
  local out=$1

  mkdir -p beds
  curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | gzip -cd - | cut -f 6-8 > beds/rMask.bed
  curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz | gzip -cd - | cut -f 2-4 > beds/trf.bed
  curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz | gzip -cd - | awk '{print $2"\t"$3"\t"$4}' > beds/gaps.bed
  cat beds/*.bed | sort -k1,1 -k2,2n > all.bed
	bedtools maskfasta -fi $fasta -bed all.bed -fo $out
  rm -f all.bed
}

kmer_locations() {
  while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
    echo "$src_dir/makeIntervalsBED_v2.pl $chrm $len 36 5" | submit -s $chrm
  done < $chrm_info
}

extract_kmers() {
  local fasta=$1
  while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
    bed="${chrm}_k36_step5_intervals.bed"
    out="${chrm}.reads.fa"
    echo "fastaFromBed -fi $fasta -bed $bed -fo $out" | submit -s $chrm
  done < $chrm_info

}

map_kmers() {
  local fasta=$1
  while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    reads="${chrm}.reads.fa"
    out="${chrm}.sam"
    echo "mrsfast --search $fasta --seq $reads -o $out -e 2" | submit -s $chrm -m 10G
  done < $chrm_info
}

count_hits() {
  local fasta=$1
  while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    kmer_file="${chrm}.reads.fa"
    alignments="${chrm}.sam"
    cmd="$src_dir/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl \
      $alignments $kmer_file $chrm > counts.$chrm"
    echo "$cmd" | submit -s count.$chrm -m 10G
  done < $chrm_info
}

mask() {
  local fasta=$1
  local out=$2
  echo "cat counts* > all.counts; maskFastaFromBed -fi $fasta -bed all.counts -fo $out" | \
  submit -s mask -m 10G
}

# Main
#######
chrm_info="chrm_info.bed"
fa_anno="mask_rtf_simple_repeat_gaps.fa"
fa_without_over="anno_kmer_masked.fa"

# Download the chrm sizes for the reference genome
gen_chrm_info

# Download annotations and apply them
mask_rtf_simple_gaps $fa_anno; index $fa_anno

kmer_locations

extract_kmers $fa_anno

map_kmers $fa_anno

count_hits

mask $fa_anno $fa_without_over
