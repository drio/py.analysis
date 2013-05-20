#!/bin/bash
set -e

TEST=0
if [ $TEST -eq 1 ];then
  locations="small.locations.bed"
  bams="."
else
  locations="locations.bed"
  bams="/stornext/snfs6/rogers/drio_scratch/baboon.diversity/papio_papio_analysis/attempt2/bams"
  if [ ! -f $locations ];then
    gzip -cd merged.vcf.gz | grep -v "#" | awk '{print $1"\t"$2-1"\t"$2}' > $locations
  fi
fi

for b in $bams/*.bam;do
  id=`echo $b | ruby -ne 'puts $_.match(/(\d+)\.bam/)[1]'`
  samtools depth -q1 -b $locations $b | sed 's/\t/_/' | sort -k1,1 -S10G > $id.depth.txt &
  #echo "$cmd" | submit -s one -l nodes=1:ppn=2,mem=8Gb
done
wait
join *.depth.txt > coverage.txt

#echo "cat small.locations.bed | coverageBed -abam $bams/papio_28547.bam -b stdin > 28547.coverage.txt " | submit -s one -l nodes=1:ppn=2,mem=8Gb
#echo "cat small.locations.bed | coverageBed -abam $bams/papio_30388.bam -b stdin > 30388.coverage.txt"  | submit -s two -l nodes=1:ppn=2,mem=8Gb
