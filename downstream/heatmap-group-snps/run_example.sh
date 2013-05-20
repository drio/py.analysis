#!/bin/bash
set -e
min_qual=$1
input_file=$2

[ ".$min_qual" == "." ] && min_qual=200
if [ ".$input_file" == "." ];then
  echo "`basename $0` <min_qual> <input_vcf_file>"
  exit 1
fi

echo "min_qual = $min_qual" >&2

rm -f *.png
(
grep "#" $input_file ; \
awk -v m=$min_qual '{if ($6>m) print;}' \
  $input_file
) | tee _tmp | wc -l

cat _tmp | ./heat-grp-snps.py test_data/groups.oconnor.tsv test_data/mhc_haplotype.tsv
#mv heat.png heat-pc-groupping.png

#cat _tmp | ./heat-grp-snps.py test_data/mhc_haplotype.tsv
#mv heat.png heat-mhc-groupping.png

rm -f _tmp
