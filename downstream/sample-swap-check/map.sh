#!/bin/bash
#
# vim: set ts=2 sw=2 noet paste:
#
ref=$1
sample=$2
[ ! -f $ref ] && echo "I can't find $ref" && exit 1
[ ".$sample" == "." ] && echo "ERROR[`basename`]: Need sample name as second parameter" && exit 1
# TODO: check for data in stdin

cat - > $fq.$sample
bwa samse $ref <(bwa aln -t4 $ref $fq.$sample) $fq.$sample | samtools view -Shb /dev/stdin > alignments.$sample.bam
rm -f $fq.$sample
samtools flagstat alignments.$sample.bam > stats.$sample.txt
