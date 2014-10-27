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
src_dir=$(dirname $(readlink -f $0))

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
#_cmd="fastqc -o . $input_bam"
#echo -e "$_cmd\tfastqc\t-"
#echo



_tmp="/space1/tmp"
_out="/dev/stdout"

# NOTE: We want to keep the original strand of the read, that
# is why we use fastqtosam instead of samtools.
_one="java -Xmx14G -jar $PICARD/SamToFastq.jar \
        TMP_DIR=$_tmp INPUT=$input_bam \
        FASTQ=one.fq SECOND_END_FASTQ=two.fq"
#_one="samtools view - "

_three="$src_dir/subreads.py | split -d -l 20000000 - ${sample_id}.fq. "
_four="gzip *fq*; rm -f one.txt two.txt"
_cmd="$_one ; cat one.txt two.txt | $_three ; $_four"
echo  "$_cmd" | submit -s subreads -m 8G -c 2
