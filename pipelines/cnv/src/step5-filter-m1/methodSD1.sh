#!/bin/bash
# vim: ts=2 expandtab:

# METHOD M1 FOR CALLING DUPLICATIONS
#
# Selects 1-Kbps windows with copy number exceeding the sample-specific cutoff
# (but lower than 100 copies), merging adjacent windows into single regions
#
# Selects regions with at least five windows and larger than 10 Kbps
#
# Retains duplications whose at least 85% of their size does not overlap with repeats
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -t 0 ]
then
  echo "Ups! I need mrcanavar output in stdin."
  exit 1
fi

for b in R sortBed mergeBed intersectBed groupBy subtractBed
do
  `which $b &>/dev/null` || error "$b not in path."
done

# Prepare a bed file with the regions considered as repeats, gaps, etc...
#
BED_RMASK=./rMask.bed # TODO: ??????????????
# curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | \
#   gzip -c - | sed 's/chr//g' | grep -e _repeat |  awk '{print $6"\t"$7"\t"$8}' > simple_repeats.bed
BED_TRF=./simple_repeats.bed
[ ! -f $BED_RMASK ] && echo "Need $BED_RMASK" && exit 1
[ ! -f $BED_TRF ] && echo "Need $BED_TRF" && exit 1
BED_REPEATS=repeats.rMask.trf.bed # merges all previous files in this one

# curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz |\
# awk '{print $2"\t"$3"\t"$4}' | sed 's/chr//g' > gaps.bed
BED_GAPS=gaps.bed
[ ! -f $BED_GAPS ] && echo "Need $BED_GAPS" && exit 1

# Computations start here
# Merge repeats in one file
#
cat $BED_RMASK $BED_TRF | sortBed -i stdin | mergeBed -i stdin > $BED_REPEATS

# Compute cuttoffs
#
R CMD BATCH $SRC_DIR/copyNumberDistribution.R
gain=`cat sample.specific.cutoffs.txt |  tail -1 | cut -f6`
echo $sample "gain cutoff: " $gain

# Compute the stretches of CNV
#
dups=./final.calls.dups.m1.bed
# In stdin we have the calls from canavar
cat - | \
  sed 1,2d | \
  grep -v 'chrY\|chrX\|chrUn\|chrM'| \
  awk -v gain=$gain '{if (gain <= $5 && $5 < 100) print}' | \
  cut -f 1,2,3,5 | \
  mergeBed -i stdin -n | \
  awk '{if ($4>=5 && $3-$2>10000) print}' | \
  intersectBed -wo -a stdin -b $BED_REPEATS | \
  groupBy -i stdin -g 1,2,3,4 -c 8 -o sum |  \
  awk '{OFS="\t"; if ($5 < 0.85*($3-$2)) print $1,$2,$3}' | \
  cut -f 1,2,3 > $dups

dups_woGaps=./final.calls.dups.m1.woGaps.bed
subtractBed -a $dups -b $BED_GAPS > $dups_woGaps
