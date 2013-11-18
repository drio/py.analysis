#!/bin/bash

mkdir -p merge
cd merge

INPUT=""
PWD=`pwd`
for b in $PWD/../sampe/*.bam
do
    INPUT="$INPUT INPUT=$b "
done

java -Xmx4g -jar $PICARD/MergeSamFiles.jar \
    $INPUT \
    SORT_ORDER=coordinate \
    USE_THREADING=True \
    OUTPUT=merged.sorted.bam
