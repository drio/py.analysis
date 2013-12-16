#!/bin/bash

TMP_DIR=$1
MEM=$2
ID=$3
URL=$4

mkdir -p dups
cd dups

INPUT="../merge/merged.sorted.bam"
java -Xmx${MEM} -jar $PICARD/MarkDuplicates.jar \
    TMP_DIR=$TMP_DIR \
    INPUT=$INPUT \
    METRICS_FILE=metrics.txt \
    OUTPUT=${ID}.merged.sorted.dups.bam

signal.py $URL $ID dups

rm -f ../merge/*.bam
