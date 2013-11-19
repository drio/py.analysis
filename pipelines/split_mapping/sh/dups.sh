#!/bin/bash

mkdir -p dups
cd dups

INPUT=`ls ../merge/*.bam`
java -Xmx4g -jar $PICARD/MarkDuplicates.jar \
    INPUT=$INPUT \
    METRICS_FILE=metrics.txt \
    OUTPUT=merged.sorted.dups.bam
