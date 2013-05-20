#!/bin/bash
#
#
reads=$1
ref=$2
n_threads=$3
probes_name=$4

usage() {
  local msg=$1
  echo $msg
  echo "$0 <reads.fastq> <path_to_ref_fasta> <n_threads> <probes_name>"
  exit 1
}

[ ".$reads" == "." ] && usage "Need reads fastq."
[ ".$ref" == "." ] && usage "Need fasta ref file."
[ ".$n_threads" == "." ] && usage "Need number of threads"
[ ".$probes_name" == "." ] && usage "Need probes name"

bwa aln -t$n_threads $ref $reads > $probes_name.sai
bwa samse $ref $probes_name.sai $reads > $probes_name.sam
rm -f $probes_name.sai
