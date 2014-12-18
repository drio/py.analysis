#!/bin/bash
# vim:set ts=2 sw=2 et paste foldmethod=indent:
#

src_dir=$(dirname $(readlink -f $0))
bin="$src_dir/../../bin"


error() {
  local msg=$1
  [ ".$msg" != "." ] && echo "ERROR: $msg"
  exit 1
}


usage() {
  local msg=$1
  [ ".$msg" != "." ] && echo  "ERROR: $msg"
  echo "Usage: $(basename $0) <sample_id> <sample_bam> <fasta_raw> <ref_kmer_masked> <canavar_conf_bin>"
  [ ".$msg" != "." ] && exit 1
}


check_run() {
    local file=$1
    local cmd=$2
    local name=$3
    local cpus=$4
    local ram=$5

    [ ".$cpus" == "." ] && cpus=1
    [ ".$ram" == "." ] && ram=5G

    if $(ls $file >/dev/null 2>/dev/null);then
        echo "$file already there, skipping" >&2
    else
      if [ $CNV_MODE == "cluster" ];then
        echo ">> $name $cpus $ram"
        (submit -f ./deps.txt -s $name -c$cpus -m$ram "$cmd") | tee ${name}.submit | bash > ./tmp.txt
        sync; sleep 1
        mv ./tmp.txt ./deps.txt
      else
        $cmd
      fi
    fi
}


: ${CNV_MODE:?"Please, set env var CNV_MODE to single or cluster"} # single or cluster
if [ $CNV_MODE != "single" ] && [ $CNV_MODE != "cluster" ];then
  echo "CNV_MODE can only by single or cluster." 2>&1
  exit 1
fi


id=$1
input_bam=$2
fasta_raw=$3
ref_kmer_masked=$4
conf=$5
[ ".$id" == "." ] && usage "I need the sample id"
[ ! -f "$input_bam" ] && usage "I need an input bam"
[ ! -f "$fasta_raw" ] && usage "Fasta raw file not found."
[ ! -f "$ref_kmer_masked" ] && usage "Ref kmer masked not found."
[ ! -f "$conf" ] && usage "Cavavar conf bin not found."

rm -f ./deps.txt; touch ./deps.txt

cmd="$bin/bam2fastq -o ${id}#.fq $input_bam" # (#) will become _1 or _2
check_run "${id}_*.fq" "$cmd" "tofasta.$id"

cmd="$bin/bwa mem -M -t4 $fasta_raw ${id}_1.fq ${id}_2.fq | \
  $bin/samtools view -Sbh - | \
  $bin/samtools sort -@4 -O bam -T $(mktemp -p /space1/tmp) - > ${id}.bam"
echo $cmd > ./bwa_raw.sh; chmod 755 ./bwa_raw.sh
check_run ${id}.bam "./bwa_raw.sh" "bwa.$id"


cmd="java -Xmx14G -jar $bin/MergeSamFiles.jar \
  TMP_DIR=/space1/tmp \
  SORT_ORDER=coordinate \
  ASSUME_SORTED=true \
  USE_THREADING=true \
  VALIDATION_STRINGENCY=SILENT \
  INPUT=${id}.bam \
  OUTPUT=merged.sorted.bam"
check_run merged.sorted.bam "$cmd" "merge.$id" 2 16G

cmd="java -Xmx14G -jar $bin/MarkDuplicates.jar
  INPUT=merged.sorted.bam
  TMP_DIR=/space1/tmp
  METRICS_FILE=metrics.txt
  VALIDATION_STRINGENCY=SILENT
  ASSUME_SORTED=true
  CREATE_INDEX=true
  REMOVE_DUPLICATES=True
  COMPRESSION_LEVEL=9
  OUTPUT=${id}.merged.sorted.dups.bam"
check_run ${id}.merged.sorted.dups.bam "$cmd" "rdups.$id" 2 16G


# Split and generate SUBREADS
cmd="$bin/bam2fastq -o no_dups#.fq ${id}.merged.sorted.dups.bam"
cmd="$cmd; cat no_dups*.fq | $bin/cleanHeaderSpaces.py | \
  $bin/kmermaker-q_MODBEL_3 -q -i stdin -k 36 -s 36 -f 10 | \
  $bin/fastqbreak -n 15000000 -o ./${id}.fq."
echo "$cmd" > subreads.sh; chmod 755 subreads.sh
check_run "${id}.fq.*" "./subreads.sh" "split.$id"


# Map subreads against kmer masked genome (no pads)
i=1
for f in ./*.fq.*
do
  _out="${i}.sam"
  cmd="mrfast --search $ref_kmer_masked --seq $f -o ${i}.sam --outcomp -e 2"
  check_run "*.sam.gz" "$cmd" "mrfast.$id"
  i=$[$i+1]
done

# Run canavar on alignments against the pad version of the ref
cmd="mrcanavar --read --gz -conf $conf -samdir . -depth ${id}.depth "
check_run ${id}.depth "$cmd" "depth.$id"

cmd="mrcanavar --call -conf $conf -depth ${id}.depth -o ${id}.output"
check_run ${id}.output.copynumber.bed "$cmd" "calls.$id"
