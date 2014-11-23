#!/bin/bash
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

error() {
  local msg=$1
  echo $msg
  exit 1
}

ref_fasta=$1
sam_dir=$2
seed=$3

# Some sanity checks
#
[ ".$ref_fasta" == "." ] && error "Need path to refence."
[ ".$sam_dir" == "." ] && error "Need sam dir."
[ ".$seed" == "." ] && error "Need seed string (perhaps sample id?)"

[ ! -d "$sam_dir" ] && error "Sam dir not found."
for b in mrcanavar
do
  `which $b &>/dev/null` || error "$b not in path."
done

# Compute the genomic windows
touch empty
conf="conf.$seed.bin"

# TODO: We only need to do this once... we can reuse it for other samples.
[ ! -f $conf ] && \
    mrcanavar --prep --gz -fasta $ref_fasta -gaps empty -conf $conf

# Find RD (and GC content on regions for all the split alignments)
[ ! -f $seed.depth ] && \
    mrcanavar --read --gz -conf $conf -samdir $sam_dir -depth $seed.depth

# perform GC correction and call duplicated/deleted regions
[ ! -f $seed.output ] && \
    mrcanavar --call -conf $conf -depth $seed.depth -o $seed.output
