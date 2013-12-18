#!/bin/bash

set -e

TMP_DIR=$1
MEM=$2
UNMAPPED_BAM=$3
REF=$4
URL=$5
ID=$6

mkdir -p merge
cd merge

ln -s $REF
samtools view -H $UNMAPPED_BAM | grep -P "^@SQ" > ./`basename $REF`.dict

INPUT=""
PWD=`pwd`
for b in $PWD/../sampe/*.bam
do
    #INPUT="$INPUT ALIGNED_BAM=$b "
    INPUT="$INPUT INPUT=$b "
done

java -Xmx${MEM} -jar $PICARD/MergeSamFiles.jar \
    $INPUT \
    TMP_DIR=$TMP_DIR \
    SORT_ORDER=coordinate \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=LENIENT \
    OUTPUT=merged.sorted.bam

signal.py $URL $ID merge

rm -f ../sampe/*.bam

touch done.txt

