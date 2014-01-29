#!/bin/bash
#
# CNV pipeline: Second step
# INPUT: bam file and sample id
# Asses the quality of the input data (reads) and generate
# 36bp subreads from the main reads. Also, trim start/end
# of reads:
#
# 1. run fastqc
# 2. remove duplicates
# 3. generate subreads (default: 10-45 and 46-81)
#
# We could be less conservative in the trimming if the quality of
# start/end of the reads is good, but, since coverage is not
# a problem, we will trim by default.
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

error() {
  local msg=$1
  echo $msg
  exit 1
}

input_bam=$1
sample_id=$2

# Some sanity checks
#
[ ".$input_bam" == "." ] && error "Need path to input bam."
[ ".$sample_id" == "." ] && error "Need sample_id"
for b in fastqc samtools
do
  `which $b &>/dev/null` || error "$b not in path."
done
[ ! -f "$PICARD/MarkDuplicates.jar" ] && error "Picard is not in \$PICARD?"

# Fastqc
#
_cmd="fastqc -o . $input_bam"
echo -e "$_cmd\tfastqc\t-"
echo

# remove dups + generate subreads
#
_tmp="/space1/tmp"
_out="/dev/stdout"
_one="java -Xmx14G -jar $PICARD/MarkDuplicates.jar REMOVE_DUPLICATES=True TMP_DIR=$_tmp INPUT=$input_bam OUTPUT=$_out METRICS_FILE=/dev/null "
_two="samtools view - "
_three="$SRC_DIR/subreads.py | gzip - > ${sample_id}.fq.gz"
_cmd="$_one | $_two | $_three"
echo -e "$_cmd\trm_dups_gen_reads\t-"
