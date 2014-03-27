#!/bin/bash

# METHOD M1 FOR CALLING DUPLICATIONS
#
# Selects 1-Kbps windows with copy number exceeding the sample-specific cutoff
# (but lower than 100 copies), merging adjacent windows into single regions
#
# Selects regions with at least five windows and larger than 10 Kbps
#
# Retains duplications whose at least 85% of their size does not overlap with repeats
#

# Samples
samples="sample1 sample2 sample3 sample4 sample5 sample6"

BEDTOOLS_DIR="XXXXXXXXXXXXXX"

# Define paths
CURRENT_PATH=`pwd`
PROJECT_DIR=/CNV/
FINAL_CALLS=$PROJECT_DIR/Final_calls_M1
mkdir -p $FINAL_CALLS
CUTOFFS=$PROJECT_DIR/scripts/sample.specific.cutoffs.txt

#BED file with repeats from RepeatMasker and Tandem Repeat Finder
BED_GAPS=gaps.bed
BED_RMASK=rMask.bed
BED_TRF=trf.bed
BED_REPEATS=repeats.rMask.trf.bed

cat $BED_RMASK $BED_TRF | $BEDTOOLS_DIR/sortBed -i stdin | $BEDTOOLS_DIR/mergeBed -i stdin > $BED_REPEATS

for sample in $samples; do
  # Get sample-specific cutoffs
  gain=`cat $CUTOFFS | grep $sample | cut -f6`
  echo $sample "gain cutoff: " $gain

  dups=$FINAL_CALLS/$sample.final.calls.dups.m1.bed
  infile=$PROJECT_DIR/mrCNVR/$sample/calls/$sample.calls.copynumber.bed
  sed 1,2d $infile | \
    grep -v 'chrY\|chrX\|chrUn\|chrM'| \
    awk -v gain=$gain '{if (gain <= $5 && $5 < 100) print}' | \
    cut -f 1,2,3,5 | $BEDTOOLS_DIR/mergeBed -i stdin -n | \
    awk '{if ($4>=5 && $3-$2>10000) print}' | \
    $BEDTOOLS_DIR/intersectBed -wo -a stdin -b $BED_REPEATS | \
    $BEDTOOLS_DIR/groupBy -i stdin -g 1,2,3,4 -c 8 -o sum |  \
    awk '{OFS="\t"; if ($5 < 0.85*($3-$2)) print $1,$2,$3}' | \
    cut -f 1,2,3 > $dups

  dups_woGaps=$FINAL_CALLS/$sample.final.calls.dups.m1.woGaps.bed
  $BEDTOOLS_DIR/subtractBed -a $dups -b $BED_GAPS > $dups_woGaps
done

