#!/bin/bash
set -e

usage()
{
  local msg=$1
  local tool_name=$( basename ${BASH_SOURCE[0]} )
  [ ".$msg" != "." ] && echo "Ups!: $msg"
  echo "Usage: $tool_name <vcf_file.gz> <wgs|wes> <path_to_id_to_colony.csv>"
  exit 1
}

vcf=$1
[ ! -f "$vcf" ] && usage "No vcf file."

seq_type=$2
if [ "$seq_type" != "wgs" ] && [ "$seq_type" != "wes" ];then
  usage "No type."
fi
csv=$3
[ ! -f "$csv" ] && usage "No csv file."

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

COLONIES="California Oregon Wisconsin Yerkes"
FCS="NON_SYNONYMOUS_CODING START_GAINED START_LOST STOP_LOST STOP_GAINED UTR_5_PRIME raw"

# FIXME: dirty hack
for colony in $COLONIES "New England"
do
  for fc in $FCS
  do
    cmd="gzip -cd $vcf | $SRC_DIR/snps_by_colony_per_fc.sh $seq_type $csv \"$colony\" $fc > \"$colony.$fc.vcf\""
    echo $cmd | submit -s $(echo "$colony.$fc" | sed 's/\s/_/g') -m1G
  done
done

