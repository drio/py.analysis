#!/bin/bash
#
#1. find read depth for the three window schemas (Only necessary once per each genome)
#    $ mrcanavar --prep --gz -fasta ../step1/genome/rhemac2.all_chrms.fa -gaps empty -conf out.conf
#2. find read depth (+perform GC correction) on windows from #1
#    $ mrcanavar --read -conf conf.bin -samdir sam_files_are_here/ -depth mysample.depth
#3. compute absolute copy number and make the CNV calls
#    $ mrcanavar --call -conf conf.bin -depth mysample.depth -o mysample
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

# Map kmers against genome
#################################################################################
echo -e "touch empty; mrcanavar --prep --gz -fasta $ref_fasta -gaps empty -conf conf.bin\twindows\t-"
echo -e "mrcanavar --read --gz -conf conf.bin -samdir $sam_dir -depth sample.depth\trdepth\twindows"
echo -e "mrcanavar --call -conf conf.bin -depth sample.depth -o output\tcall_cnvs\trdepth"
echo



