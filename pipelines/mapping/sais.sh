#!/bin/bash

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

help() {
  local msg=$1
  [ ".$msg" != "." ] && echo "UPS!: $msg"
  echo "Usage: `basename $0` <input.bam> <ref_fasta> <read1_or_two> <output_seed> <list_of_deps> <sample_id> [n_threads]"
  exit 1
}

bam=$1
[ ".$bam" == "." ]         && help "Need input bam."
[ ! -f $bam ]              && help "Cannot find input bam."
ref_fasta=$2
[ ".$ref_fasta" == "." ]   && help "Need fasta file"
[ ! -f $ref_fasta ]        && help "Cannot find input bam."
one_or_two=$3
[ ".$one_or_two" == "." ]  && help "Tell me if I should compute read 1 or 2"
output_seed=$4
[ ".$output_seed" == "." ] && help "Need output seed"
deps=$5
[ ".$deps" == "." ]        && help "Need list of deps. Use - for none."
sample_id=$6
[ ".$sample_id" == "." ]   && help "Need sample id."
n_threads=$7
[ ".$n_threads" == "." ]   && n_threads=1

job_id=sai.$one_or_two.$sample_id.$RANDOM
output=$output_seed.$one_or_two.sai
cmd="bwa aln -t$n_threads -b$one_or_two $ref_fasta $bam > $output"
echo -e "$job_id\t5G\t$n_threads\t$deps\t$cmd\t$output"
