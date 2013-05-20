#!/bin/bash

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SRC_DIR=$ROOT_DIR/../src
REF="/stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa"

tsv="crv_wgs.tsv"

#chinese="36013|35082|35084|35086"
chinese=`cat $tsv | awk -F$'\t' '{print $2"\t"$7}' | grep -v DN | grep -i Chine | ruby -ane 'BEGIN{@a=[]}; @a << $F[0]; END{print @a.join("|")}'`
chinese_bams=`find -L ./bams -name "*.bam" | grep -P "$chinese"`
indian_bams=`find -L ./bams -name "*.bam" | grep -P -v "$chinese"`

all_bams=`find -L ./bams -name "*.bam"`
mkdir -p logs

for pf in `ls ./probes/*.fixed`
do
  seed_pf=`basename $pf`
  pfq=${pf}.readified.fq
  i=0
  INF=2000
  compute_only=$INF

  echo "## align; $pf"
  cmda="$SRC_DIR/align-probes.sh $pfq $REF 1 $seed_pf"
  echo $cmda | submit -s a.$seed_pf -m5G -c2
  #echo $cmda

  for b in $all_bams
  do
    [ $i -eq $compute_only ] && break
    i=$[$i+1]

    id=`echo $b | ruby -ne 'puts $_.match(/(\d+)/)[1]'`

    echo "## eg; $id"
    cmde="$SRC_DIR/run-eg.sh $b $pf $id $seed_pf"
    echo $cmde | submit -s e.$seed_pf.$id.$i -m4G -c2
    #echo $cmde

  done
  echo ""
done
