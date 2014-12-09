#!/bin/bash
# vim:set ts=2 sw=2 et paste foldmethod=indent:
#
# In case you want to remove:
# rm -rf anno_kmer_masked.fa *.nohit  .RData  kmerCountsWithMrsFast_more20placements.txt  moab_logs/ deps/  *.bed chr*.fa counts.chr* *.sam all*.counts
#
# We need:
: << EOF
  wget http://hgdownload.cse.ucsc.edu/goldenPath/rheMac2/bigZips/chromFa.tar.gz
  tar zxf chromFa.tar.gz
  cat *.fa > raw_fasta.fa
EOF

pipe="bash" # "" to not submit or "bash" to submit to cluster
th_over="20"
src_dir=$(dirname $(readlink -f $0))


error() {
  local msg=$1
  [ ".$msg" != "." ] && echo "ERROR: $msg"
  exit 1
}


usage() {
  local msg=$1
  [ ".$msg" != "." ] && echo  "ERROR: $msg"
  echo "Usage: $(basename $0) <raw_fasta.fa> <ucsc_genome>"
  [ ".$msg" != "." ] && exit 1
}


# Params checks
raw_fasta=$1
ucsc_genome=$2
[ ".$raw_fasta" == "." ] && usage "I need raw fasta name"
[ ".$ucsc_genome" == "." ] && usage "Need ucsc genome code"


# default filename outputs. Don't change.
fasta_anno_masked="fasta_anno_masked.fa"
fa_without_over="anno_kmer_masked.fa"
chrm_info="chrm_info.bed"


# STEP1: download chrm info for genome
########
if [ ! -f $chrm_info.bed ];then
  mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
      -e "select chrom, size from ${ucsc_genome}.chromInfo" | \
      grep chr | awk '{print $1"\t"$2"\t"$2}' | grep -v size > $chrm_info
else
  echo "$chrm_info already there."
fi

[ ! -f $chrm_info ] && error "I can't find $chrm_info" && exit 1
rm -rf deps; mkdir deps
for i in `seq 1 10`;do
  touch deps/${i}.txt
done


# STEP2: Download anotations and apply them
###################################
if [ ! -f ./$fasta_anno_masked ];then
  cmd="$src_dir/mask_rtf_simple_gaps.sh $raw_fasta $fasta_anno_masked $ucsc_genome"
  submit -e -s "gen_masked_ref" $cmd | $pipe >> deps/1.txt
else
  echo "$fasta_anno_masked already there."
fi


# STEP3: Generate read coordinates in genome
######################################
while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
		cmd="$src_dir/makeIntervalsBED_v2.pl $chrm $len 36 5"
		cat deps/1.txt | submit -f - -e -s "locations.$chrm" $cmd | $pipe >> deps/2.txt
done < $chrm_info


# STEP4: extract reads from genome
######################################
while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
    bed="${chrm}_k36_step5_intervals.bed"
    out="${chrm}.reads.fa"
    cmd="fastaFromBed -fi $fasta_anno_masked -bed $bed -fo $out"
    cat deps/2.txt | submit -e -f - -s "extract.$chrm" "$cmd" | $pipe >> deps/3.txt
done < $chrm_info


# STEP5: map genome reads to genome
####################################
while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    reads="${chrm}.reads.fa"
    out="${chrm}" # mrsfast already adds the extension .sam
    # FIXME
    # Hack to avoid failing dependency.
    # There is one read file with all N's and mrsfast exit != 0
    mrsfast="$src_dir/../../bin/mrsfast"
    cmd="$mrsfast --search $fasta_anno_masked --seq $reads -o $out -e 2; exit 0"
    cat deps/3.txt | submit -e -f - -m 16G -s "map_kmers.$chrm" "$cmd" | $pipe >> deps/4.txt
done < $chrm_info


# STEP6: Count how many alignments we have per genome location
####################################
while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    kmer_file="${chrm}.reads.fa"
    alignments="${chrm}.sam"
		# FIXME
		# Hack to avoid failing dependency.
		# There is one read file with all N's and mrsfast exits with value different to 0
    cmd="$src_dir/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl \
			$alignments $kmer_file $chrm > counts.$chrm; exit 0"
    cat deps/4.txt | submit -e -m 16G -f - -s "count.$chrm" "$cmd" | $pipe >> deps/5.txt
done < $chrm_info


# STEP7: merge counts
########################
cmd="$src_dir/process_counts.sh $th_over"
cat deps/5.txt | submit -e -f - -s "working_on_counts" $cmd | $pipe >> deps/6.txt

# STEP8: Mask overrepresented regions in the annotated fasta reference
##############################
cmd="maskFastaFromBed -fi $fasta_anno_masked -bed all.$th_over.counts -fo $fa_without_over"
cat deps/6.txt | submit -e -f - -s "mask_final" "$cmd" | $pipe
