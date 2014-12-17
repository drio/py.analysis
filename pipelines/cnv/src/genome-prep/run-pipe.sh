#!/bin/bash
# vim:set ts=2 sw=2 et paste foldmethod=indent:
#
# In case you want to remove:
# rm -rf moab_logs/ partial_output/ deps/ mrcanavar.bins  gzip_versions/*.gz empty  for_* fasta_anno_masked.fa* pads.bed  anno_kmer_masked.fa all.* anno_kmer_masked.fa.pad beds/ chrm_info.bed counts* *.sam *reads.fa *intervals* *.nohit
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
bin="$src_dir/../../bin"

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
    [ -f ${chrm}_k36_step5_intervals.bed ] && continue
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
    [ -f $out ] && continue
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
    [ -f ${out}.sam ] && continue
    cmd="$bin/mrsfast --search $fasta_anno_masked --seq $reads -o $out -e 2; exit 0"
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
    [ -f counts.$chrm ] && continue
    cmd="$src_dir/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl \
			$alignments $kmer_file $chrm > counts.$chrm; exit 0"
    cat deps/4.txt | submit -e -m 16G -f - -s "count.$chrm" "$cmd" | $pipe >> deps/5.txt
done < $chrm_info


# STEP7: merge counts
########################
if [ ! -f all.${th_over}.counts ];then
  cmd="$src_dir/process_counts.sh $th_over"
  cat deps/5.txt | submit -e -f - -s "working_on_counts" $cmd | $pipe >> deps/6.txt
fi


# STEP8: Mask overrepresented regions in the annotated fasta reference
##############################
if [ ! -f $fa_without_over ];then
  cmd="maskFastaFromBed -fi $fasta_anno_masked -bed all.${th_over}.counts -fo $fa_without_over"
  cat deps/6.txt | submit -e -f - -s "mask_final" "$cmd" | $pipe >> deps/7.txt
fi


# STEP9: extend the gaps regions with N's (padding)
if [ ! -f ${fa_without_over}.pad ];then
  cmd="$src_dir/pad_reference.sh beds/rMask.bed beds/trf.bed beds/gaps.bed ./all.${th_over}.counts $fa_without_over $chrm_info ${fa_without_over}.pad"
  cat deps/7.txt | submit -e -f - -s "padding_ref" "$cmd" | $pipe >> deps/8.txt
fi


# STEP10: compute the windows (bin file) with mrcanavar
if [ ! -f ./mrcanavar.bins ];then
  touch ./empty
  cmd="$bin/mrcanavar --prep --gz -fasta ${fa_without_over}.pad -gaps empty -conf mrcanavar.bins"
  cat deps/8.txt | submit -e -f - -s "bins" "$cmd" | $pipe
fi


# STEP11: create mrfast indexing for anno_kmer_masked reference
if [ ! -d ./for_mrfast ];then
  cmd="mkdir ./for_mrfast; cd ./for_mrfast; ln -s ../$fa_without_over; $bin/mrfast --index ./$fa_without_over; samtools faidx .$fa_without_over"
  cat deps/7.txt | submit -e -f - -s "index_mrfast" "$cmd" | $pipe
fi


# STEP12: create bwa indices for raw_reference
if [ ! -d ./for_bwa ];then
  cmd="mkdir ./for_bwa; cd ./for_bwa; ln -s ../$raw_fasta; $bin/bwa index ./$raw_fasta"
  cat deps/7.txt | submit -e -f - -s "bwa_index" "$cmd" | $pipe
fi


# STEP14: Perform sanity checks on output

