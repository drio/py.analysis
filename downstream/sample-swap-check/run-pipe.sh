#!/bin/bash
#
# vim: set ts=2 sw=2 noet paste:
#
src_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ref=$1
sample=$2
n_reads=$3

help() {
	local msg=$1
	echo "ERROR: " $msg
	echo "Usage: `basename $0` <ref.fa> <sample_name> <n_reads_to_sample>"
	[ ".$msg" != "." ] && exit 1
}

[ ! -f $ref ] && help "I can't find $ref"
[ ".$sample" == "." ] && help "Need sample name as second parameter"
if test -t 0
then
	help "Need fastq stream in stdin"
fi

cat /dev/stdin |\
	$src_dir/down-sample.sh $n_reads |\
	$src_dir/map.sh $ref $sample
