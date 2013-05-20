#!/bin/bash

# gzip -cd merged.anno.rdp.vcf.gz| grep -v "#" | grep -v INDEL | ruby -ne 'puts $_.match(/EFF=([\w_]+)/)[1]' | sort  | uniq -c | sort -k1,1n
FCS="SYNONYMOUS_START NON_SYNONYMOUS_START SYNONYMOUS_STOP STOP_LOST START_LOST STOP_GAINED START_GAINED UTR_5_PRIME UTR_3_PRIME NON_SYNONYMOUS_CODING"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
  local msg=$1
  [ ".$msg" != "." ] && echo "ERROR: $msg"
  cat <<EOF
Usage: cat file.vcf | tool <FUNC_CONS>"
for fc in $FCS;do submit -s \$fc "gzip -cd merged.anno.vcf.gz | $DIR/allele.freq.plots.sh \$fc"; done
EOF
  exit 1
}

[ -t 0 ] && usage "No data in stdin"
[ $# -ne 1 ] && usage "Wrong number of args"
FC=$1
[ ".$FC" == "." ] && usage "Need functional consequence."

# AN: Total number of alleles in called genotypes
data=data.$FC.txt
if [ ! -f $data ];then
  cat - | \
    grep -v "#" | \
    grep -v INDEL | \
    #grep PASS | \
    awk '{if ($6>50) print}' | \
    grep -P "$FC" | \
    ruby -ne 'puts $_.match(/AN=([\w_]+)/)[1]' | \
    sort -n -T /tmp/ | uniq -c | sort -k1,1n \
    > $data
fi

png=af.$FC.png
(
cat <<EOF
set t png size 800,200
set title "Variant Allele frequencies in $FC substitutions"
set style line 1 lt 1 lw 8
set xlabel "Number of events"
set ylabel "log(Allele Counter)"
set logscale y
plot "-" using 2:1 notitle with imp ls 1
EOF
cat $data
) | gnuplot > $png
ls -lachd $png
