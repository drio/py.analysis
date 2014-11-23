#!/bin/bash
#
usage() {
    local msg=$1
    [ ".$msg" != "." ] && echo "$msg"
    echo "`basename $0` <input_bam> <n_reads_per_split>"
    [ ".$msg" != "." ] && exit 1 || exit 0
}

[ $# -ne 2 ] && usage "Wrong number of parameters."
input_bam=$1
n_reads=$2
[ ! -f $input_bam ] && usage "I can't find input_bam."
[ "$n_reads" == "."] && usage "Need # of reads per split."

split -d -l $n_reads - split.fq.
