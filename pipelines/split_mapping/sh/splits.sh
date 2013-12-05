#!/bin/bash
set -e

INPUT_BAM=$1
NR_X_SPLIT=$2
MEM=$3
TMP_DIR=$4

[ `uname` == "Darwin" ] && SPLIT="gsplit" || SPLIT="split"

# We have to sort by query name otherwise we
# may end up breaking two pairs and causing
# trouble down the pipeline. Specifically,
# this will happen:
# http://sourceforge.net/mailarchive/forum.php?thread_name=50993D76.8060400%40broadinstitute.org&forum_name=samtools-help
# Where, you are missing one of the ends of read pair.
#
mkdir -p splits
cd splits
samtools view -H $INPUT_BAM > header.sam
samtools view -H $INPUT_BAM | grep -P "^@RG"> rg.sam
samtools view $INPUT_BAM | \
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

ls *.bam | wc -l > ./done.txt
