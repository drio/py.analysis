#!/bin/bash
#
set -e

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage()
{
  echo "Usage: snp.freq.sh [wgs|wes]"
  exit 1
}

seq_type=$1
[ "$seq_type" != "wgs" ] && [ "$seq_type" != "wes" ] && usage
[ $seq_type == "wes" ] && gen_size=34000000 || gen_size=2900000000
single_freq="0"
n_samples=0
for i in `find . -name "*.vcf.anno.gz" | xargs -i echo {}`
do
  if [ "$i" != "./merged.vcf.anno.gz" ];then # skip the merged one
    num=`gzip -cd $i | $SRC_DIR/single.filter.sh $seq_type | wc -l`
    id=`echo $i | ruby -ne 'puts $_.match(/\/(\d+)\./)[1]'`
    echo -ne "$id : "
    bc=`echo "scale=10;($num/$gen_size)*1000" | bc`
    single_freq="$single_freq $bc"
    n_samples=$[n_samples+1]
    echo "scale=10;($single_freq)/$n_samples" | sed 's/ /\+/g' |bc
  fi
done

if [ $n_samples -eq 0 ];then
  echo "I couldn't process any vcf file. Are you in the right dir?"
else
  echo "scale=10;($single_freq)/$n_samples" | sed 's/ /\+/g' |bc
fi
