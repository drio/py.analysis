#!/bin/bash
#
[ ! -f "../common.sh" ] && echo "I can't find ../common.sh" && exit 1
[ ! -f "./chrm_info.bed" ] && echo "I can't find ./chrm_info.bed" && exit 1
source ../common.sh

N_READS_SPLIT=5000000 # /2

# Build the script we will use to process split reads and only
# output the ones that are not properly(1) map by bwa.
# Do this to send as less reads as possible to mrfast (it is slow).
(
  cat <<EOF
#!/bin/bash
#
source ../common.sh

bwa aln $fasta \$1 |  \
  bwa samse $fasta - \$1 |  \
  grep -v "*$" |  \
  awk '{if (\$5 == 0) print ">"\$1"\n"\$10;}'
EOF
) > bwa.sh; chmod 755 ./bwa.sh

split() {
  local chrm=$1
  local len=$2
  local dep=$3
  cmd="$step1_sh intervals  -c $chrm -l $len | $step1_sh extract -f $fasta -b - | split -d -l $N_READS_SPLIT - split.${chrm}."
  echo $cmd | submit -s $chrm.split
}

bwa_reduce() {
  for s in split.*
  do
    cmd="./bwa.sh $s > tmp.$s; mv tmp.$s $s"
    echo $cmd | submit -s bwa.$s -m6g
  done
}

map_counts() {
  local chrm=$1
  local len=$2
  local dep=$3
  for s in split.${chrm}.*
  do
    # FIXME: very broken, stdout does not stream as expected
    #cmd="$step1_sh map -f $fasta -r bwa.$s | $step1_sh count -m bwa.$s -a - -c $chrm > counts.$s"
    cmd="mrsfast --search $fasta --seq $s -o /dev/stdout -e 2 | `dirname $step1_sh`/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl /dev/stdin $s $chrm > counts.$s"
    echo $cmd | submit -s map.$s -m 8g
  done
}

mask() {
  cat counts* | $step1_sh mask -f $fasta -t - > masked.fa
  cat counts* | xvfb-run R CMD BATCH --no-restore --no-save '--args stdin' \
    /stornext/snfs7/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/step1-masking/makeHistogramCounts.R
}

# All chrms
all_chrms() {
  local step=$1
  local dep=$2
  while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
    $step $chrm $len $dep
  done < chrm_info.bed
}


index() {
	echo "mrsfast --index $fasta --ws 12"
}

# Main
#######
echo "Start here..."
# Mask rtf, simple, gaps

#index
#all_chrms "split"
#bwa_reduce
#all_chrms "map_counts" "null"
#all_chrms "mask"
