#!/bin/bash

TMP_DIR=$1
MEM=$2
UNMAPPED_BAM=$3
REF=$4

mkdir -p merge
cd merge

ln -s $REF
samtools view -H $UNMAPPED_BAM | grep -P "^@SQ" > ./`basename $REF`.dict

INPUT=""
PWD=`pwd`
for b in $PWD/../sampe/*.bam
do
    INPUT="$INPUT ALIGNED_BAM=$b "
done

java -Xmx${MEM}g -jar $PICARD/MergeBamAlignment.jar \
    $INPUT \
    UNMAPPED_BAM=$UNMAPPED_BAM \
    REFERENCE_SEQUENCE=./`basename $REF` \
    PAIRED_RUN=true \
    TMP_DIR=$TMP_DIR \
    SORT_ORDER=coordinate \
    OUTPUT=merged.sorted.bam
