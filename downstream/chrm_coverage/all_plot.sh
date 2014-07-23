#!/bin/bash

for cov_f in *.gz
do
  for chrm in `cat genome.bed  | awk '{print $1}' | grep -v random`
  do
    sample=`echo $cov_f | sed 's/.bam.bed.gz//g'`
    submit -s ${sample}_${chrm}_png -m8G "gzip -cd $cov_f | grep -P \"^${chrm}\t\" | ./plot.py ${sample}_${chrm} > ${sample}.${chrm}.png"
  done
done

