#!/bin/bash

FASTA=$1
SAI1=$2
SAI2=$3
SPLIT_NUMBER=$4
BAM=$5
MEM=$6

mkdir -p sampe
cd sampe

bwa sampe $FASTA $SAI1 $SAI2 $BAM $BAM | \
    java -Xmx${MEM}g -jar $PICARD/SortSam.jar \
    SORT_ORDER=coordinate INPUT=/dev/stdin OUTPUT=$SPLIT_NUMBER.bam
