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

# http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_MarkDuplicates_.28or_ValidateSamFile.29_produces_the_error_.22Value_was_put_into_PairInfoMap_more_than_once..22__What_should_I_do.3F
java -Xmx${MEM}g -jar $PICARD/MergeBamAlignment.jar \
    $INPUT \
    UNMAPPED_BAM=$UNMAPPED_BAM \
    REFERENCE_SEQUENCE=./`basename $REF` \
    PAIRED_RUN=true \
    TMP_DIR=$TMP_DIR \
    SORT_ORDER=coordinate \
    OUTPUT=merged.sorted.bam
