#!/bin/bash
#
# vim: set ts=2 sw=2 noet paste:
#
src_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ref=$1
sample=$2
n_reads=$3
input_fq=$4

help() {
	local msg=$1
	echo "ERROR: " $msg
	echo "Usage: `basename $0` <ref.fa> <sample_name> <n_reads_to_sample> <input_fq>"
	[ ".$msg" != "." ] && exit 1
}

[ ! -f $ref ] && help "I can't find $ref"
[ ! -f $input_fq ] && help "I can't find $input_fq"
[ ".$sample" == "." ] && help "Need sample name as second parameter"
# [ -t 0 ] && help "Need fastq stream in stdin"

# cat - | 
# $src_dir/down-sample.sh $n_reads  | 
# $src_dir/map.sh $ref $sample

r_string=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
seed=$sample.$r_string


echo "Sampling ..."
cat $input_fq | \
	awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | \
	shuf  | \
	head -$n_reads | \
	sed 's/\t\t/\n/g' > $seed.fq


echo "Mapping ..."
bwa aln -t4 $ref $seed.fq > $seed.sai
bwa samse $ref $seed.sai $seed.fq > $seed.sam 
samtools view -Shb $seed.sam > alignments.$seed.bam
samtools flagstat alignments.$seed.bam > stats.$sample.$n_reads.$r_string.txt

echo "Cleaning up ..."
rm -f $seed.fq alignments.$seed.bam $seed.sai $seed.sam
