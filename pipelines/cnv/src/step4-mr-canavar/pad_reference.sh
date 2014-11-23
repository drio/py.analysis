#!/bin/bash

error() {
  local msg=$1
  [ ".$msg" != "." ] && echo "ERROR: $msg"
  echo "`basename $0` <r_master> <tr_finder> <gaps_bed> <kmer_bed> <input_ref> <output>"
  exit 1
}

r_master=$1
tr_finder=$2
gaps_bed=$3
kmer_bed=$4
input_ref=$5
output=$6

[ ! -f $r_master ] && error "I can't find repeat masker bed"
[ ! -f $tr_finder ] && error "I can't find tr_finder masker bed"
[ ! -f $gaps_bed ] && error "I can't find gaps_bed masker bed"
[ ! -f $kmer_bed ] && error "I can't find kmer_bed masker bed"
[ ! -f $input_ref ] && error "I can't find input_ref"
[ ".$output" == "."  ] && error "I need output file name"

cat $r_masker $tr_finder $gaps_bed $kmer_bed | \
  awk '{print $1"\t"$2"\t"$3}' | \
  sortBed | \
  mergeBed | \
  awk '{
    start=$2; 
    if (start-36>0) { start=start-36; } 
    end=$3+36; 
    print $1"\t"start"\t"end
    }' |\
  bedtools maskfasta -fi $input_ref -bed stdin -fo $output
