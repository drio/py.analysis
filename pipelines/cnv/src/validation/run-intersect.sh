#!/bin/bash
# vim: ts=2 expandtab ft=sh:
#
set -e

usage() {
  local _exit=$1
  cat <<EOF
Usage:
  $ tool <cutoff_raw_data_min> <cutoff_raw_data_max> <cutoff_truth> <cutoff_truth_max>
EOF
  exit $_exit
}

error() {
  local msg=$1
  echo "ERROR: $msg"
  usage 1
}

cutoff_raw_data_min=$1
cutoff_raw_data_max=$2
cutoff_truth_min=$3
cutoff_truth_max=$4

[ ".$cutoff_raw_data_min" == "." ] && error "Need cutoff min param for raw set."
[ ".$cutoff_raw_data_max" == "." ] && error "Need cutoff max param for raw set."
[ ".$cutoff_truth_min" == "." ] && error "Need cutoff min param for truth set."
[ ".$cutoff_truth_min" == "." ] && error "Need cutoff max param for truth set."

[ ! -d ./truth ] && error "No truth dir found."
[ ! -d ./filter ] && error "No filter dir found."

run_first_part=true
if $run_first_part
then
  cd truth

#  echo "Preparing merged.bed" 1>&2
#  cat bak/*.bed | sed 's/chr//' | grep -vP "^X|^Y" | sed 's/10+/10/' | \
#    awk -v c=$cutoff_truth '{if ($4>=c) print}' | \
#    cut -f 1,2,3 | sort -k1,1n -k2,2n -T=/tmp -S4G | \
#    mergeBed -i stdin > merged.bed

  for i in bak/*bed
  do
    o=`basename $i | sed 's/Homo_sapiens-//g' | sed -s 's/_1000bp_simple_per_1kb.HM.bedgraph.bed//g'`
    echo $o >&2
    cat $i | sed 's/chr//' | grep -vP '^X|^Y' | sed 's/10+/10/' | \
      awk -v a=$cutoff_truth_min -v c=$cutoff_truth_max '{if ($4>=a && $4<=c) print}' | sort -k1,1n -k2,2n > $o &
  done
  wait

  cd ..
fi

cd filter
# Our Raw calls
#
echo "Preparing raw calls" >&2
cut -f1,2,3,5 output.copynumber.bed | grep -vP "^X|^Y" | sed '1,2d' | \
  awk -v a=$cutoff_raw_data_min -v c=$cutoff_raw_data_max '{if ($4>=a && $4<=c) print;}' | sort -k1,1n -k2,2n > raw_calls.bed

# Perform the intersections
#
our_calls_file="raw_calls.bed"
bp_raw=`awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}' $our_calls_file`
bp_raw_si=`echo $bp_raw | numfmt --to=si`

for i in ../truth/*HGD*
do
  _bn=`basename $i`

  bp_if=`awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}' $i`
  bp_if_si=`echo $bp_if | numfmt --to=si`

  # if -> i
  bp_int=`intersectBed -u -wa -a $our_calls_file -b $i | awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}'`
  _p1=`echo "scale=2; ($bp_int/$bp_raw)*100" | bc`
  echo "$our_calls_file,$_bn,$bp_raw_si,$bp_if_si,$_p1"

  # i -> if
  bp_int=`intersectBed -u -wa -a $i -b $our_calls_file | awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}'`
  _p1=`echo "scale=2; ($bp_int/$bp_if)*100" | bc`
  echo "$_bn,$our_calls_file,$bp_if_si,$bp_raw_si,$_p1"
done
cd ..
