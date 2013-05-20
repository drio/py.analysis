#!/bin/bash
set -e

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
N_CONCURRENT=$1
LOCATIONS=$2
COVERAGE=$3
INPUT_VCF=$4

# Compute concurrently the coverage per each bam
MAPQ=10
find . -maxdepth 1 -name "*.bam"| xargs -t -I {} -n 1 -P $N_CONCURRENT sh -c \
  "samtools depth -q$MAPQ -b $LOCATIONS {} | gzip -c > ./{}.depth.txt.gz"

# each line will contain the coverage at locus per each sample
$SRC_DIR/my_gnu_join.py *.depth.txt.gz | gzip -c > $COVERAGE

# Method2: broken

#state=first
#_tmp=$TMP_SORT/_tmp
#_cov=$TMP_SORT/_cov
#for f in `gzip -cd $INPUT_VCF | grep "^#CHROM" | ruby -ane 'puts $F[9..-1].join("*.depth.txt ") + "*.depth.txt"'`
#do
#  if [ $state == 'first' ];then
#    prev=$f
#    state='second'
#  elif [ $state == 'second' ];then
#    state='normal'
#    echo "join $prev $f >$_cov" >&2
#    join $prev $f > $_cov
#  else
#    echo "join $_cov $f > $_tmp; mv $_tmp $_cov" >&2
#    join $_cov $f > $_tmp
#    mv $_tmp $_cov
#  fi
#done
#gzip -c $_cov > $COVERAGE
#rm -f $_cov

#echo "cat small.locations.bed | coverageBed -abam $bams/papio_28547.bam -b stdin > 28547.coverage.txt " | submit -s one -l nodes=1:ppn=2,mem=8Gb
#echo "cat small.locations.bed | coverageBed -abam $bams/papio_30388.bam -b stdin > 30388.coverage.txt"  | submit -s two -l nodes=1:ppn=2,mem=8Gb
