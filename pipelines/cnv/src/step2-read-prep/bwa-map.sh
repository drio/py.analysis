#!/bin/bash
#
usage() {
    local msg=$1
    [ ".$msg" != "." ] && echo "$msg"
    echo "`basename $0` <genome.fa> <input_bam> > output.bam"
    [ ".$msg" != "." ] && exit 1 || exit 0
}

[ $# -ne 2 ] && usage "Wrong number of parameters."
fa=$1
input_bam=$2
[ ! -f $fa ] && usage "I can't find genome."
[ ! -f $input_bam ] && usage "I can't find input_bam."

# Notice that the bwa aln run in parallel
# Plan accordingly.
bwa sampe $fa \
  <(bwa aln -b -1 $fa $input_bam) \
  <(bwa aln -b -2 $fa $input_bam) \
  $input_bam $input_bam |\
samtools view -Shb /dev/stdin


