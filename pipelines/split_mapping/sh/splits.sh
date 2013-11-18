#!/bin/bash

INPUT_BAM=$1
NR_X_SPLIT=$2

[ `uname` == "Darwin" ] && SPLIT="gsplit" || SPLIT="split"

mkdir -p splits
cd splits
samtools view -H $INPUT_BAM > header.sam
samtools view $INPUT_BAM | \
    $SPLIT -d -l $NR_X_SPLIT - split.
for s in split*
do
    (cat header.sam ; cat $s) | samtools view -Sb - > $s.bam
done

ls split.* | grep -v bam | xargs rm -f

ls *.bam | wc -l > ./done.txt
