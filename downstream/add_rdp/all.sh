#!/bin/bash
#
#set -e

error() {
  (
  echo "ERROR: $1"
  cat <<-EOF
Usage: $ $TOOL <bams dir> <input_vcf_file> <num of concurrent processes> <test 0|1>
EOF
  ) >&2
  exit 1
}

show_env_vars() {
  (
  cat <<-EOF
---------------
bams=$BAMS_DIR
input=$INPUT_VCF
cores=$N_CONCURRENT
testing=$TEST
----------------
EOF
  ) >&2
}

proccess_args() {
  BAMS_DIR=$1
  INPUT_VCF=$2
  N_CONCURRENT=$3
  TEST=$4

  [ ! -d $BAMS_DIR ]         && error "Bams dir not found"
  [ ! -f $INPUT_VCF ]        && error "Cannot find input vcf."
  OUTPUT_VCF="$INPUT_VCF.rdp.gz"
  [ "$N_CONCURRENT" == "." ] && error "Need # of cores."
  [ "$TEST" == "1" ] && HEAD="head -10000" || HEAD="cat -"
  SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  show_env_vars
}

dowork() {
  if [ ! -f $LOCATIONS ]
  then
    echo "`date`>> Creating locations HEAD='$HEAD'" >&2
    gzip -cd $INPUT_VCF | grep -v "#" | awk '{print $1"\t"$2-1"\t"$2}' | $HEAD > $LOCATIONS
  fi

  # Find coverage per each sample
  echo "`date`>> Joining Coverage n_cores=$N_CONCURRENT" >&2
  MAPQ=10
  find $BAMS_DIR -maxdepth 1 -name "*.bam" | xargs -t -I {} -n 1 -P $N_CONCURRENT sh -c \
    "[ ! -f {}.depth.txt.gz ] && samtools depth -q$MAPQ -b $LOCATIONS {} |\
      $SRC_DIR/add_missing_locations.py $LOCATIONS | gzip -c > ./{}.depth.txt.gz || echo 'skipping.'"

  # Join coverage at locus for all samples
  # NOTE: the order in which we feed the files to my_gnu_join.py is important and has to match
  # the order of the samples in the vcf header
  depth_files=`gzip -cd $INPUT_VCF | grep "^#CHROM" | ruby -ane 'puts $F[9..-1].join("*.depth.txt.gz ") + "*.depth.txt.gz "'`
  [ ! -f $COVERAGE ] && $SRC_DIR/my_gnu_join.py $depth_files | gzip -c > $COVERAGE || \
    echo 'skipping.'

  echo "`date`>> Adding RDP" >&2
  [ ! -f $OUTPUT_VCF ] && gzip -cd $COVERAGE | \
    $SRC_DIR/add_rdp.py $INPUT_VCF | gzip -c > $OUTPUT_VCF || echo 'skipping.'
}

########
# Main
########
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TOOL="add_rdp.sh"
COVERAGE="coverage.txt.gz"
LOCATIONS="locations.bed"

# Defaults
INPUT_VCF=""
BAMS_DIR=""
TEST=0
N_CONCURRENT=`cat /proc/cpuinfo  | grep proc | wc -l`

proccess_args $1 $2 $3 $4
dowork

#./final_answer.sh
#ssh drio@is04607.com rm -f public_html/tmp/papi*
#scp papio.vcf.anno.rdp.q200.rdp10.mapq10.gz drio@is04607.com:public_html/tmp/
#echo "http://is04607.com/tmp/papio.vcf.anno.rdp.q200.rdp10.mapq10.gz"
