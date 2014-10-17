#!/bin/bash
# vim: ts=2 expandtab ft=sh:
#
set -e

usage() {
  local _exit=$1
  cat <<EOF
Usage:
  $ tool <cutoff_raw_data_min> <cutoff_raw_data_max> <cutoff_truth> <cutoff_truth_max> <min_win_size> <truth_bed> <calls_bed>

Example:
  $ ./src/run-intersect.sh 3.5 20 3 10 1000 inputs/french.bed inputs/hg18_drio_canavar.bed
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
min_win_size=$5
truth_file=$6
calls_file=$7

[ ".$cutoff_raw_data_min" == "." ] && error "Need cutoff min param for raw set."
[ ".$cutoff_raw_data_max" == "." ] && error "Need cutoff max param for raw set."
[ ".$cutoff_truth_min" == "." ] && error "Need cutoff min param for truth set."
[ ".$cutoff_truth_min" == "." ] && error "Need cutoff max param for truth set."
[ ".$truth_file" == "." ] && error "Need truth file."
[ ".$calls_file" == "." ] && error "Need calls file."

RAND=`strings /dev/urandom | grep -o '[[:alnum:]]' | head -n 30 | tr -d '\n'; echo`

echo "Preparing truth ($truth_file)" >&2
# chr1    0       1231    2
p_truth="truth.$RAND"
sed 's/chr//' $truth_file | grep -vP '^X|^Y|^GL|random' | sed 's/10+/10/' | \
  awk -v a=$cutoff_truth_min -v c=$cutoff_truth_max '{if ($4>=a && $4<=c) print}' | \
  awk -v ws=$min_win_size '{if ($3-$2>ws) print}' | \
  sort -k1,1 -k2,2n > $p_truth

echo "Preparing raw calls ($calls_file)" >&2
#CHROM  START   END     GC%     COPYNUMBER
#
#chr1    0       2679    0.577000        21.883331
our_calls_file="raw_calls.bed.$RAND"
n_cols=`awk '{print NF}' $calls_file | sort -nu | tail -n 1`
if [ $n_cols -eq 5 ]; then
  echo "<5> columns detected." >&2
  sed 's/chr//' $calls_file | cut -f1,2,3,5 | grep -vP "^X|^Y|^GL|random" | sed '1,2d' | \
    awk -v a=$cutoff_raw_data_min -v c=$cutoff_raw_data_max '{if ($4>=a && $4<=c) print;}' | \
    awk -v ws=$min_win_size '{if ($3-$2>ws) print}' | \
    sort -k1,1 -k2,2n > $our_calls_file
elif [ $n_cols -eq 3 ]; then
  echo "<$n_cols> columns detected." >&2
  sed 's/chr//' $calls_file | grep -vP "^X|^Y|^GL|random" |  \
    awk -v ws=$min_win_size '{if ($3-$2>ws) print}' | \
    sort -k1,1 -k2,2n > $our_calls_file
elif [ $n_cols -eq 4 ]; then
  echo "<4> columns detected." >&2
  sed 's/chr//' $calls_file | grep -vP "^X|^Y|^GL|random" | sed 's/10+/10/' | \
    awk -v a=$cutoff_raw_data_min -v c=$cutoff_raw_data_max '{if ($4>=a && $4<=c) print;}' | \
    awk -v ws=$min_win_size '{if ($3-$2>ws) print}' | \
    sort -k1,1 -k2,2n > $our_calls_file
else
  echo "Incorrect number of columns in input calls. Bailing out."
  exit 1
fi

# Perform the intersections
#
num_raw_calls=`cat $our_calls_file | wc -l`
bp_raw=`awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}' $our_calls_file`
bp_raw_si=`echo $bp_raw | numfmt --to=si`

# Compute metrics
i=$p_truth
_bn=`basename $i`
num_truth_calls=`cat $i | wc -l`

bp_truth_calls=`awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}' $i`
bp_if_si=`echo $bp_truth_calls | numfmt --to=si`

# intersect calls vs truth
num_int=`intersectBed -u -wa -a $our_calls_file -b $i | wc -l`
bp_int=`intersectBed -u -wa -a $our_calls_file -b $i | awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}'`
_p1=`echo "scale=2; ($bp_int/$bp_raw)*100" | bc`
_n1=`echo "scale=2; ($num_int/$num_raw_calls)*100" | bc`

# Fraction of our calls that are miss-calls
_mc=`echo "scale=2; ($num_raw_calls-$num_int)/$num_raw_calls" | bc`
_mc_bp=`echo "scale=2; ($bp_raw-$bp_int)/$bp_raw" | bc`

# intersect truth vs calls
num_int=`intersectBed -u -wa -a $i -b $our_calls_file | wc -l`
bp_int=`intersectBed -u -wa -a $i -b $our_calls_file | awk 'BEGIN{a=0}; {a = ($3-$2) + a}; END{print a}'`
_p2=`echo "scale=2; ($bp_int/$bp_truth_calls)*100" | bc`
_n2=`echo "scale=2; ($num_int/$num_truth_calls)*100" | bc`

# Fraction of our calls that are in the truth set
_sensitivity=`echo "scale=2; $num_int/$num_truth_calls" | bc`
_sensitivity_bp=`echo "scale=2; $bp_int/$bp_truth_calls" | bc`


#echo "input1,input2,bp1,bp2,bp_overlap_1-2,bp_overlap2-1,num_1-2,num2-1,cutoff_1_min,cutoff_1_max,cutoff_2_min,cutoff_2_max,min_win_size,sensitivity,miss_rate,sensitivity_bp,miss_rate_bp" | sed 's/,/\t/g' >&2
#echo "$our_calls_file,$_bn,$bp_raw_si,$bp_if_si,$_p1,$_p2,$_n1,$_n2,$cutoff_raw_data_min,$cutoff_raw_data_max,$cutoff_truth_min,$cutoff_truth_max,$min_win_size,$_sensitivity,$_miss_rate,$_sensitivity_bp,$_miss_rate_bp" | sed 's/,/\t/g'
echo "nc1,nc2,bp1,bp2,bp_overlap_1-2,bp_overlap2-1,num_1-2,num2-1,cutoff_1_min,cutoff_1_max,cutoff_2_min,cutoff_2_max,min_win_size,\
  sensitivity,sensitivity_bp,mc,mc_bp" | sed 's/,/\t/g' >&2
echo "$num_raw_calls,$num_truth_calls,$bp_raw_si,$bp_if_si,$_p1,$_p2,$_n1,$_n2,$cutoff_raw_data_min,$cutoff_raw_data_max,$cutoff_truth_min,$cutoff_truth_max,$min_win_size,\
  $_sensitivity,$_sensitivity_bp,$_mc,$_mc_bp" | sed 's/,/\t/g'

rm -f $p_truth $our_calls_file
