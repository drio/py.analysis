#!/bin/bash
#
src_dir="/stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/cnv"
bams_dir="`pwd`/bams"
ref_bam="$bams_dir/35087.bam"
list_chrms=`samtools view -H $ref_bam | grep SN | awk -F: '{print $2}' | awk '{print $1}' | grep -v random`
x_coverage=30
read_length=100

window_size=$1
if [ ".$window_size" == "." ];then
  echo "Need window size"
  exit 1
fi

output_dir="log2ratios/$window_size"
num_reads_per_window=`echo "$x_coverage * ($window_size/$read_length)" | bc`

mkdir -p $output_dir
ref_id=`basename $ref_bam | sed 's/.bam//g'`
for b in $bams_dir/*.bam; do
  id=`basename $b | sed 's/.bam//g'`
  for chrm in $list_chrms; do
    seed=${ref_id}_vs_${id}_${chrm}
    o_file=$output_dir/${seed}.txt.gz
    cmd="$src_dir/binner.py $num_reads_per_window $ref_bam $chrm | $src_dir/compute_log_ratio.py $b | gzip -c > $o_file"
    echo "$cmd" | submit -s "$seed" -c3 -m4
  done
done
