#!/bin/bash

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

help() {
  local msg=$1
  [ ".$msg" != "." ] && echo "UPS!: $msg"
  cat << EOF
Usage: `basename $0` <options>
  -b <input.bam>
  -f <ref_fasta>
  -o <output_seed>
  -r <ram_picard>
  -m <tmp_dir>
  -s <sample_id>
  -t [n_threads_sampe]
  -d // to disable smith waterman in sampe

Example:
  $ $HOME/dev/py.analysis/pipelines/mapping/jobs_for_mapping.sh -b /Users/drio/dev/bam.examples/phix.bam -f /Users/drio/dev/genomes/phix.fa -o FOOOOO -r 1G -m /tmp -s 11833 -t 4 | awk -F\t '{print \$5}'  | grep -v cmd | bash
EOF
  exit 1
}

n_threads=1 # defaults
disable_sw="false"
while getopts "b:f:o:r:m:s:t:d" opt
do
  case $opt in
    b)  bam=$OPTARG ;;
    f)  ref_fasta=$OPTARG     ;;
    o)  output_seed=$OPTARG   ;;
    r)  ram=$OPTARG           ;;
    m)  tmp=$OPTARG           ;;
    s)  sample_id=$OPTARG     ;;
    d)  disable_sw="true"    ;;
    t)  n_threads=$OPTARG     ;;
    \?) help "Invalid option: -$OPTARG" ;;
    :) help "Option -$OPTARG requires an argument." ;;
  esac
done

[ ".$bam" == "." ]         && help "Need input bam."
[ ! -f $bam ]              && help "Cannot find input bam."
[ ".$ref_fasta" == "." ]   && help "Need fasta file"
[ ! -f $ref_fasta ]        && help "Cannot find ref fasta."
[ ".$output_seed" == "." ] && help "Need output seed"
[ ".$ram" == "." ]         && help "Need amount of ram in Gb. Example: 4"
[ ".$tmp" == "." ]         && help "Need tmp dir"
[ ".$sample_id" == "." ]   && help "Need sample id."

# Main
#
#r_string=$(rand_string)
jobs_file=job_list.$output_seed.tsv

$SRC_DIR/sais.sh $bam $ref_fasta 1 $output_seed - $sample_id $n_threads > $jobs_file
sai1_id=`cat $jobs_file | tail -1 | awk -F"\t" '{print $1}'`
sai1_out=`cat $jobs_file | tail -1 | awk -F"\t" '{print $6}'`
$SRC_DIR/sais.sh $bam $ref_fasta 2 $output_seed - $sample_id $n_threads >> $jobs_file
sai2_id=`cat $jobs_file | tail -1 | awk -F"\t" '{print $1}'`
sai2_out=`cat $jobs_file | tail -1 | awk -F"\t" '{print $6}'`

$SRC_DIR/sampe.sh $bam $ref_fasta $sai1_out $sai2_out $output_seed $sample_id $disable_sw $sai1_id,$sai2_id >> $jobs_file
sampe_id=`cat $jobs_file | tail -1 | awk -F"\t" '{print $1}'`
sampe_out=`cat $jobs_file | tail -1 | awk -F"\t" '{print $6}'`

$SRC_DIR/bam_coordiante_sort.sh $sampe_out $output_seed.sorted.bam $ram $tmp $sampe_id $sample_id >> $jobs_file
sort_id=`cat $jobs_file | tail -1 | awk -F"\t" '{print $1}'`
sort_out=`cat $jobs_file | tail -1 | awk -F"\t" '{print $6}'`

$SRC_DIR/mark_dups.sh $sort_out $output_seed.sorted.dups.bam $ram $tmp $sample_id $sort_id >> $jobs_file
dups_id=`cat $jobs_file | tail -1 | awk -F"\t" '{print $1}'`
dups_out=`cat $jobs_file | tail -1 | awk -F"\t" '{print $6}'`

echo -e "id\tram\tthreads\tdeps\tcmd\toutput_files"
echo -e "rm.$RANDOM\t1G\t1\t$dups_id\trm -f $sai1_out $sai2_out $sampe_out $sort_out\t-" >> $jobs_file

cat $jobs_file
rm -f $jobs_file

