#!/bin/bash

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

help() {
  local msg=$1
  [ ".$msg" != "." ] && echo "UPS!: $msg"
  echo "Usage: `basename $0` <input_bam> <output_bam> <ram> <tmp_dir> <sample_id> [list_of_deps]"
  exit 1
}

bam=$1
[ ".$bam" == "." ]         && help "Need input bam."
output=$2
[ ".$output" == "." ]      && help "Need output bam"
ram=$3
[ ".$ram" == "." ]         && help "Need amount of ram in Gb. Example: 4G"
tmp=$4
[ ".$tmp" == "." ]         && help "Need tmp dir"
deps=$5
[ ".$deps" == "." ]        && deps="-"
sample_id=$6
[ ".$sample_id" == "." ]   && help "Need sample id."

job_id=2bam.$sample_id.$RANDOM
cmd="java -Xmx$ram -jar $PICARD/SortSam.jar SORT_ORDER=coordinate TMP_DIR=$tmp INPUT=$bam OUTPUT=$output VALIDATION_STRINGENCY=LENIENT"
echo -e "$job_id\t$ram\t1\t$deps\t$cmd\t$output"



