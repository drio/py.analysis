#!/bin/bash
#
# samtools view -h ../bams/18277.bam | head -1000000 | samtools view -Sb - > 18277.bam
set -e

#BAM=$1
#[ ".$BAM" == "." ] && usage

BAM=./18277.bam
#if [ ! -f $BAM ];then
  echo ">> Copying bam from ardmore"
  AMOUNT=500000
  #ssh ardmore "samtools view -h /stornext/snfs6/rogers/drio_scratch/rhesus/crv.v2/wgs/human-array-validation/bams/$BAM |\
  samtools view -h /stornext/snfs6/rogers/drio_scratch/rhesus/crv.v2/wgs/human-array-validation/bams/$BAM |\
  head -$AMOUNT |\
  samtools view -bS - > $BAM
#fi

#REF=/Users/drio/dev/genomes/hsap_36.1_hg18.normal.80chr.fa
#REF=/Users/drio/dev/genomes/rhemac2.indian_macaque_no_phix.fixed.fa
REF=/stornext/snfs6/rogers/drio_scratch/genomes/hsap_36.1_hg18.normal.80chr.fa

# -O3       : Gap open penalty [11]
# -E3       : Gap extension penalty [4]
# -R1000000 : Proceed with suboptimal alignments if there are no more than INT equally best hits.
# -i1       : -i INT   Disallow an indel within INT bp towards the ends [5]
echo ">> bwa aln"
bwa aln -O3 -E3 -R1000000 -i1 -t4 -b1 $REF $BAM > sai
echo ">> bwa samse"
bwa samse $REF sai $BAM > out.sam
cat out.sam | java -jar $PICARD/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT=out.bam
#rm -f sai
