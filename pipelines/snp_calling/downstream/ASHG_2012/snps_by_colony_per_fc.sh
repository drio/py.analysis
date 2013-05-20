#!/bin/bash
#

usage()
{
  local msg=$1
  local tool_name=$( basename ${BASH_SOURCE[0]} )
  [ ".$msg" != "." ] && echo "Ups!: $msg"
  echo "Usage: cat some.vcf | $tool_name <wgs|wes> <ids_to_colonies.csv> <colony_name> <func_cons>"
  echo "Use func_cons = raw if don't want to filter at all."
  exit 1
}

seq_type=$1
csv=$2
colony=$3
func_cons=$4

if [ "$seq_type" != "wgs" ] && [ "$seq_type" != "wes" ];then
  usage "No type."
fi
[ ! -f "$csv" ] && usage "No csv file."
[ ".$colony" == "." ] && usage "No colony name."
if [ "$func_cons" == "raw" ];then
  this_only="cat -"
else
  this_only="grep EFF=$func_cons"
fi
echo "this_only = $this_only" >&2

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
source $SRC_DIR/common.sh

arr=("$csv" "$colony")
cat - |\
  $SRC_DIR/colony_filter.py "${arr[@]}" |\
  $SRC_DIR/population.filter.sh $seq_type | $this_only
