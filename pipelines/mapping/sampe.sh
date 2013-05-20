#!/bin/bash

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

help() {
  local msg=$1
  [ ".$msg" != "." ] && echo "UPS!: $msg"
  echo "Usage: `basename $0` <input_bam> <ref_fasta> <sai1> <sai2> <output_seed> <sample_id> <disable_sw> [list_of_deps]"
  exit 1
}

bam=$1
[ ".$bam" == "." ]         && help "Need input bam."
[ ! -f $bam ]              && help "Cannot find input bam."
ref_fasta=$2
[ ".$ref_fasta" == "." ]   && help "Need fasta file"
[ ! -f $ref_fasta ]        && help "Cannot find input bam."
sai1=$3
[ ".$sai1" == "." ]        && help "Need sai1 file"
sai2=$4
[ ".$sai2" == "." ]        && help "Need sai2 file"
output_seed=$5
[ ".$output_seed" == "." ] && help "Need output seed"
sample_id=$6
[ ".$sample_id" == "." ]   && help "Need sample id."
disable_sw=$7
[ ".$disable_sw" == "." ]  && help "Need to tell me if you want to disable sw (true) or not (false)"

deps=$8
[ ".$deps" == "." ]        && deps="-"

job_id=sampe.$sample_id.$RANDOM
output=$output_seed.sam
[ `echo $disable_sw | tr '[:lower:]' '[:upper:]'` == "TRUE" ] && sw_o="-s"
cmd="bwa sampe $sw_o $ref_fasta $sai1 $sai2 $bam $bam > $output"
echo -e "$job_id\t5G\t1\t$deps\t$cmd\t$output"


