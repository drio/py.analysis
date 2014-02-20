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

# Some sanity checks
#
[ ".$ref_fasta" == "." ] && error "Need path to refence."
[ ".$sam_dir" == "." ] && error "Need sam dir."
[ ! -d "$sam_dir" ] && error "Sam dir not found."
for b in mrcanavar
do
  `which $b &>/dev/null` || error "$b not in path."
done

# Compute the genomic windows
touch empty
echo -e "mrcanavar --prep --gz -fasta $ref_fasta -gaps empty -conf conf.bin\twindows\t-"
echo

# Find RD (and GC content on regions for all the split alignments)
depth_files=""
for sam in $sam_dir/*.sam.gz
do
    bn=`basename $sam`
    mkdir -p ./$bn
    [ ! -f ./$bn/$bn ] && ln -s `readlink -e $sam` ./$bn/$bn
    single_sam_dir=`pwd`/$bn
    echo -e "mrcanavar --read --gz -conf conf.bin -samdir $single_sam_dir -depth ${bn}.depth\trdepth.${bn}\twindows"
    depth_files="$depth_files ${bn}.depth"
done
echo

# Merge the results from the previous steps
echo -e "mrcanavar --conc -conf conf.bin -concdepth $depth_files -depth merged.depth\tmergeRD\trdepth"
echo

# perform GC correction and call duplicated/deleted regions
echo -e "mrcanavar --call -conf conf.bin -depth merged.depth -o output\tcall_cnvs\tmergeRD"
echo
