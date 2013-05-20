#!/bin/bash
#
set -e

usage()
{
  echo "Usage: single.filter_vcf [wgs|wes]"
  exit 1
}

seq_type=$1
if [ "$seq_type" != "wgs" ] && [ "$seq_type" != "wes" ];then
  usage
fi

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/common.sh

# Population snp vcf
# Filter by annotation Func cons (wes), quality, min number of samples
cat - | \
  grep -v "#" | grep -v INDEL | \
  awk -v qual=$min_snp_qual '{if ($6>=qual) print}' |\
  ruby -ne 'ns=$_.match(/AN=(\d+);/)[1].to_i/2; puts $_ if ns > ENV["min_num_samples"].to_i'
