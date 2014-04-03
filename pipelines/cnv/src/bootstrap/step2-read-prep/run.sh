#!/bin/bash
#
#
source ../common.sh

#samtools view -h $bam | head -100000 |
#  java -Xmx4G -jar /stornext/snfs6/rogers/drio_scratch/bb/local/picard/SortSam.jar SORT_ORDER=coordinate INPUT=/dev/stdin OUTPUT=foo.bam

if [ ".$1" == "." ]
then
  $step2_sh $bam $id | awk -F"\t" '{print $1}' | grep fastq | submit -s fq.$id
  $step2_sh $bam $id | awk -F"\t" '{print $1}' | grep java | submit -s trim.$id -m16G
else
  [ ! -f foo.bam ] && \
  samtools view -h $bam | head -100000 |
   java -Xmx4G -jar /stornext/snfs6/rogers/drio_scratch/bb/local/picard/SortSam.jar SORT_ORDER=coordinate INPUT=/dev/stdin OUTPUT=foo.bam;
  /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/step2-read-prep/run-pipe.sh ./foo.bam foo | awk -F"\t" '{print $1}' | grep java
fi
