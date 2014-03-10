#!/bin/bash

bam=$1

if [ ".$bam" == "." ];then
  echo "Need input bam."
  exit 1
fi

[ ! -f ./genome.bed ] && samtools view -H $bam  | grep SQ | awk -F: '{print $2"1\t"$3}' | sed 's/LN//g' > genome.bed
/hgsc_software/BEDTools/latest/bin/bedtools genomecov -g genome.bed -ibam $bam -d | gzip -c > cov.$bam.bed.gz
