#!/bin/bash

TMP_DIR=$1
MEM=$2

mkdir -p merge
cd merge

INPUT=""
PWD=`pwd`
for b in $PWD/../sampe/*.bam
do
    INPUT="$INPUT INPUT=$b "
done

java -Xmx${MEM}g -jar $PICARD/MergeSamFiles.jar \
    $INPUT \
    TMP_DIR=$TMP_DIR \
    SORT_ORDER=coordinate \
    USE_THREADING=True \
    OUTPUT=merged.sorted.bam
