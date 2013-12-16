#!/bin/bash
set -e

INPUT_BAM=$1
NR_X_SPLIT=$2
MEM=$3
TMP_DIR=$4
URL=$5
ID=$6


[ `uname` == "Darwin" ] && SPLIT="gsplit" || SPLIT="split"

mkdir -p splits
cd splits
samtools view -H $INPUT_BAM > header.sam
samtools view -H $INPUT_BAM | grep -P "^@RG"> rg.sam

# We need to sort the input bam by queryname so we have
# the two reads of a pair adjacent to each other.
samtools view -h $INPUT_BAM | \
java -Xmx${MEM} -jar $PICARD/SortSam.jar \
    INPUT=/dev/stdin \
    TMP_DIR=$TMP_DIR \
    SORT_ORDER=queryname \
    OUTPUT=/dev/stdout | \
        samtools view - | \
        $SPLIT -d -l $NR_X_SPLIT - split.

for s in split*
do
    (cat header.sam ; cat $s) |  \
    java -Xmx${MEM} -jar $PICARD/SortSam.jar \
        INPUT=/dev/stdin \
        TMP_DIR=$TMP_DIR \
        SORT_ORDER=queryname \
        OUTPUT=$s.bam
done

ls split.* | grep -v bam | xargs rm -f

signal.py $URL $ID splits

ls *.bam | wc -l > ./done.txt
