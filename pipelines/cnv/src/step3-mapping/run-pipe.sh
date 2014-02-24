#!/bin/bash
#
# CNV pipeline: Third step
# INPUT: splited fastq file (.gz)
# Map the splited files back to the reference
#
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

error() {
  local msg=$1
  echo $msg
  exit 1
}

input_dir=$1
ref_fasta=$2

# Some sanity checks
#
[ ".$input_dir" == "." ] && error "Need path to fastq files."
[ ".$ref_fasta" == "." ] && error "Need path to refence."
for b in mrsfast
do
  `which $b &>/dev/null` || error "$b not in path."
done

# Map kmers against genome
#################################################################################
i=1
for f in $input_dir/*.gz
do
  _out="${i}.sam"
  echo -e "mrsfast --search $ref_fasta --seq $f -o $_out --seqcomp --outcomp -e 2\tmapping\t-"
  i=$[$i+1]
done
echo




