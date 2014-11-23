#!/bin/bash
#
# vim: set ts=2 sw=2 noet paste:
#
ref=$1
sample=$2

[ ! -f $ref ] && echo "I can't find $ref" && exit 1
[ ".$sample" == "." ] && echo "ERROR[`basename`]: Need sample name as second parameter" && exit 1
# TODO: check for data in stdin

r_string=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
seed=$sample.$r_string

cat /dev/stdin > $seed.fq
n_lines=$(wc -l $seed.fq | awk '{print $1}')
n_reads=$(echo "$n_lines/4" | bc)

bwa aln -t4 $ref $seed.fq > $seed.sai
bwa samse $ref $seed.sai $seed.fq > $seed.sam 
samtools view -Shb $seed.sam > alignments.$seed.bam
samtools flagstat alignments.$seed.bam > stats.$sample.$n_reads.$r_string.txt
rm -f $seed.fq alignments.$seed.bam $seed.sai $seed.sam
