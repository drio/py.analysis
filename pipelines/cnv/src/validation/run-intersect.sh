#!/bin/bash
# vim: ts=2 expandtab ft=sh:
#
set -e

usage() {
  local _exit=$1
  cat <<EOF
Usage:
  $ tool <cutoff_raw_data_min> <cutoff_raw_data_max> <cutoff_truth> <cutoff_truth_max> <truth_bed> <calls_bed>
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
truth_file=$5
calls_file=$6


[ ".$cutoff_raw_data_min" == "." ] && error "Need cutoff min param for raw set."
[ ".$cutoff_raw_data_max" == "." ] && error "Need cutoff max param for raw set."
[ ".$cutoff_truth_min" == "." ] && error "Need cutoff min param for truth set."
[ ".$cutoff_truth_min" == "." ] && error "Need cutoff max param for truth set."
[ ".$truth_file" == "." ] && error "Need truth file."
[ ".$calls_file" == "." ] && error "Need calls file."


echo "Preparing truth" >&2
p_truth="truth.$RANDOM.$RANDOM"
sed 's/chr//' $truth_file | grep -vP '^X|^Y|^GL|random' | sed 's/10+/10/' | \
  awk -v a=$cutoff_truth_min -v c=$cutoff_truth_max '{if ($4>=a && $4<=c) print}' | sort -k1,1n -k2,2n > $p_truth

echo "Preparing raw calls" >&2
n_cols=`awk '{print NF}' $calls_file | sort -nu | tail -n 1`
if [ $n_cols -eq 5 ]; then
  echo "<5> columns detected." >&2
  cut -f1,2,3,5 $calls_file | grep -vP "^X|^Y|^GL|random" | sed 's/chr//' | sed '1,2d' | \
    awk -v a=$cutoff_raw_data_min -v c=$cutoff_raw_data_max '{if ($4>=a && $4<=c) print;}' | sort -k1,1n -k2,2n > raw_calls.bed
elif [ $n_cols -eq 3 ]; then
  echo "<3> columns detected." >&2
  grep -vP "^X|^Y|^GL|random" $calls_file | sed 's/chr//' | sort -k1,1n -k2,2n > raw_calls.bed
else
  echo "Incorrect number of columns in input calls. Bailing out."
  exit 1
fi

# Perform the intersections
#
our_calls_file="raw_calls.bed"
bp_raw=`awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}' $our_calls_file`
bp_raw_si=`echo $bp_raw | numfmt --to=si`

for i in $p_truth
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

rm -f $p_truth
