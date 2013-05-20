#!/bin/bash
#
set -e

usage()
{
  echo "Usage: single.filter_vcf <wgs|wes> [drop_fcs]"
  echo "if "drop_fcs" passed; this FCs will be filtered out: $f_cons_to_remove"
  exit 1
}

seq_type=$1
if [ "$seq_type" != "wgs" ] && [ "$seq_type" != "wes" ];then
  usage
fi

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/common.sh

cat - | \
  grep -v "#" | $drop_this | grep -v INDEL | \
  awk -v qual=$min_snp_qual '{if ($6>=qual) print}' |\
  ruby -ne 'dp=$_.match(/DP=(\d+);/)[1].to_i; puts $_ if dp > ENV["min"].to_i && dp < ENV["max"].to_i'
