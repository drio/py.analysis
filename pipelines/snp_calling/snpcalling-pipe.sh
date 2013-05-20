#!/bin/bash
#
src_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

help() {
  local msg=$1
  [ ".$msg" != "." ] && echo "ERROR: $msg"
  echo "help: tool <bams_dir> <wes_or_wgs> <ref_fasta> <snp_eff_id>"
  exit 1
}

[ $# -lt 3 ] && help "Wrong number of args"
b_dir=$1
[ ".$b_dir" == "." ] && help "Need bam dir."
wes_or_wgs=$2
[ ".$wes_or_wgs" == "." ] && help "Is this a wes or wgs sample?"
ref_fasta=$3
[ ".$ref_fasta" == "." ] && help
[ ! -f "$ref_fasta" ] && help "I couldn't open the fasta file."
snp_eff_db=$4

for b in $b_dir/*.bam
do
  id=`echo $b | ruby -ne 'puts $_.match(/([\w_-]+)\.bam/)[1]'`
  if [ -f "$id.bam" ];then
    echo "# $id.bam is already there. Skipping."
  else
    submit -c2 -s $id "$src_dir/single-snpcall.sh $id $b $wes_or_wgs $ref_fasta $snp_eff_db | bash"
  fi
done
