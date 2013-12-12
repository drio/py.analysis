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
bwa sampe -r "`cat ../splits/rg.sam`" -P -A $FASTA $SAI1 $SAI2 $BAM $BAM > $TMP_DIR/${SPLIT_NUMBER}.sam

java -Xmx${MEM} -jar $PICARD/FixMateInformation.jar \
    INPUT=$TMP_DIR/${SPLIT_NUMBER}.sam OUTPUT=${SPLIT_NUMBER}.bam TMP_DIR=$TMP_DIR

rm -f $SAI1 $SAI2
rm -f $TMP_DIR/${SPLIT_NUMBER}.sam
rm -f $BAM
