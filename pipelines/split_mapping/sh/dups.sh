#!/bin/bash

set -e

TMP_DIR=$1
MEM=$2
ID=$3

mkdir -p dups
cd dups

INPUT="../merge/merged.sorted.bam"
java -Xmx${MEM} -jar $PICARD/MarkDuplicates.jar \
    TMP_DIR=$TMP_DIR \
    INPUT=$INPUT \
    METRICS_FILE=metrics.txt \
    OUTPUT=${ID}.merged.sorted.dups.bam

rm -f ../merge/*.bam
