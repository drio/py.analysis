#!/bin/bash
#
#
reads=$1
probes=$2
id=$3
probes_name=$4

usage() {
  local msg=$1
  echo $msg
  echo "$0 <reads.bam> <probes.txt> <sample id> <probes_name>"
  exit 1
}

[ ".$reads" == "." ]  && usage "Need reads (bam)."
[ ".$probes" == "." ] && usage "Need probes file."
[ ".$id" == "." ]     && usage "Need id."
[ ".$probes_name" == "." ] && usage "Need probes name"

#( java -jar $PICARD/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT=$reads FASTQ=/dev/stdout SECOND_END_FASTQ=/dev/null; \
#  java -jar $PICARD/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT=$reads FASTQ=/dev/null SECOND_END_FASTQ=/dev/stdout) |\
samtools view -F 1024 $reads |\
  ruby -ane 'puts ">#{$F[0]}\n#{$F[9]}"' |\
  egrl hits -p $probes -r - |\
  awk -F, '{print $1}' |\
  sort -T$TMPDIR -S2G -k1,1 |\
  uniq -c > $id.counts.$probes_name.txt
