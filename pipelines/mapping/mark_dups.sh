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
sample_id=$5
[ ".$sample_id" == "." ]   && help "Need sample id."

deps=$6
[ ".$deps" == "." ]        && deps="-"

job_id=dups.$sample_id.$RANDOM
metrics=$RANDOM.metrics
cmd="java -Xmx$ram -jar $PICARD/MarkDuplicates.jar INPUT=$bam OUTPUT=$output METRICS_FILE=$metrics TMP_DIR=$tmp VALIDATION_STRINGENCY=LENIENT; rm -f $metrics"
echo -e "$job_id\t$ram\t1\t$deps\t$cmd\t$output"


