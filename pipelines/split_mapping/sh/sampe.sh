#!/bin/bash

FASTA=$1
SAI1=$2
SAI2=$3
SPLIT_NUMBER=$4
BAM=$5
MEM=$6
TMP_DIR=$7

mkdir -p sampe
cd sampe

# TODO: set -A as an option
bwa sampe -r "`cat ../splits/rg.sam`" -P -A $FASTA $SAI1 $SAI2 $BAM $BAM | \
java -Xmx${MEM} -jar $PICARD/SortSam.jar \
    TMP_DIR=$TMP_DIR SORT_ORDER=coordinate \
    INPUT=/dev/stdin OUTPUT=$TMP_DIR/tmp_${SPLIT_NUMBER}.bam
java -Xmx${MEM} -jar $PICARD/FixMateInformation.jar \
    INPUT=$TMP_DIR/tmp_${SPLIT_NUMBER}.bam OUTPUT=$SPLIT_NUMBER.bam TMP_DIR=$TMP_DIR
rm -f $TMP_DIR/tmp_${SPLIT_NUMBER}.bam
