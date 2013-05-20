#!/bin/bash

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DEF_REF="$HOME/dev/bam.examples/PhiX_plus_SNP.fa"

usage() {
  cat <<EOF
`basename $0` <fastq file for the probes> <input_bam> <sample_id> <reference.fa> <probes_files> <probe_nameh>
Example:
  $ `basename $0` ./test_data/probes.fq test_data/18277.bam 18277 $DEF_REF ./test_data/affy6.fixed.txt affy
  $ `basename $0` /stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/human-array-validation/test_data/probes.fq /stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/human-array-validat.fa /stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/human-array-validation/test_data/affy6.fixed.txt affy6
EOF
}

error() {
  local what=$1
  echo "Need input parameter: $what"
  usage
  exit 1
}

FQ_PROBES=$1
BAM=$2
ID=$3
REF=$4
PROBES=$5
PROBES_NAME=$6

[ ! -f "$FQ_PROBES" ] && error "Cannot find fq probe file."
[ ! -f "$BAM" ]    && error "Cannot find bam file."
[ ! -f "$REF" ]    && error "Cannot find reference file."
[ ! -f "$PROBES" ] && error "Cannot find probe file."
[ ".$ID" == "." ]   && error "Id not provided"
[ ".$PROBES_NAME" == "." ]   && error "need probes name."

# fix probes
#cat $PROBES | ../src/fix_probe_file.sh > affy6.fixed.txt
# Readify probes
#cat $PROBES | ../src/probes2reads.py > probes.fq
# align readifyied probes to the reference genome
$SRC_DIR/align-probes.sh $FQ_PROBES $REF 4 $PROBES_NAME
# egenotype the raw reads against the regular probe file
$SRC_DIR/run-eg.sh $BAM $PROBES $ID $PROBES_NAME
