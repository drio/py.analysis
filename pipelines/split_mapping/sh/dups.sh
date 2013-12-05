#!/bin/bash

set -e

TMP_DIR=$1
MEM=$2

mkdir -p dups
cd dups

INPUT=`ls ../merge/*.bam`
java -Xmx${MEM} -jar $PICARD/MarkDuplicates.jar \
    TMP_DIR=$TMP_DIR \
    INPUT=$INPUT \
    METRICS_FILE=metrics.txt \
    OUTPUT=merged.sorted.dups.bam

#rm -f ../merge/*.bam
#rm -f ../sampe/*.bam
#rm -f ../splits/*.bam
