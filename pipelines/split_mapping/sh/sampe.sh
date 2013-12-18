#!/bin/bash

set -e

FASTA=$1
SAI1=$2
SAI2=$3
SPLIT_NUMBER=$4
BAM=$5
MEM=$6
TMP_DIR=$7
URL=$8
ID=$9



mkdir -p sampe
cd sampe

# TODO: set -A as an option
#bwa sampe -r "`cat ../splits/rg.sam`" -P -A $FASTA $SAI1 $SAI2 $BAM $BAM > $TMP_DIR/${SPLIT_NUMBER}.sam

bwa sampe -r "`cat ../splits/rg.sam`" -P -A $FASTA $SAI1 $SAI2 $BAM $BAM | \
    java -Xmx${MEM} -jar $PICARD/FixMateInformation.jar \
        VALIDATION_STRINGENCY=LENIENT \
        INPUT=/dev/stdin OUTPUT=${SPLIT_NUMBER}.bam TMP_DIR=$TMP_DIR

signal.py $URL $ID sampe

rm -f $SAI1 $SAI2
rm -f $BAM
