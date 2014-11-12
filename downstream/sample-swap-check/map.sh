#!/bin/bash
#
# vim: set ts=2 sw=2 noet paste:
#
ref=$1
sample=$2

[ ! -f $ref ] && echo "I can't find $ref" && exit 1
[ ".$sample" == "." ] && echo "ERROR[`basename`]: Need sample name as second parameter" && exit 1
# TODO: check for data in stdin

r_string=$(openssl rand -base64 8)
seed=$sample.$r_string

cat - > $seed.fq
n_lines=$(wc -l $seed.fq | awk '{print $1}')
n_reads=$(echo "$n_lines/4" | bc)

bwa samse $ref <(bwa aln -t4 $ref $seed.fq) $seed.fq | samtools view -Shb /dev/stdin > alignments.$seed.bam
samtools flagstat alignments.$seed.bam > stats.$sample.$n_reads.$r_string.txt
rm -f $seed.fq alignments.$seed.bam
