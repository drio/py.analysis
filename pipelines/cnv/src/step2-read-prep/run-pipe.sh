#!/bin/bash
#
# vim: set ts=2 sw=2 et:
#
src_dir=$(dirname $(readlink -f $0))

error() {
  local msg=$1
  echo $msg
  exit 1
}


for b in fastqc samtools bam2fastq
do
  `which $b &>/dev/null` || error "$b not in path."
done


input_bam=$1
sample_id=$2

[ ".$input_bam" == "." ] && error "Need path to input bam."
[ ".$sample_id" == "." ] && error "Need sample_id"

cp $src_dir/*.c $src_dir/clean*.py .
gcc kmermaker-q_MODBEL_3.c -o ./kmermaker-q_MODBEL_3
gcc -lm fastqbreak.c -o ./fastqbreak  &>/dev/null

[ ! -f ./kmermaker-q_MODBEL_3 ] && error "I can't compile kmermaker"
[ ! -f ./fastqbreak ] && error "I can't compile fastqbreak"

cmd="bam2fastq -o no_dups.#.fq $input_bam"

cmd="$cmd; cat no_dups*.fq | ./cleanHeaderSpaces.py | \
  ./kmermaker-q_MODBEL_3 -q -i stdin -k 36 -s 36 -f 10 | \
  ./fastqbreak -n 15000000 -o ./$sample_id.fq"

# cmd="$cmd cat no_dups*.fq | $src_dir/subreads.py | split -d -l 200000000 - ${sample_id}.fq."
cmd="$cmd; gzip $sample_id.fq*"

echo "$cmd" | submit -s subreads -m 12G -c 2
