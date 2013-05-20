#!/bin/bash

MIN_STD=0.07
MIN_NUM_ALL_PER_GROUP=10
main_dir=$HOME/dev/py.analysis

export PYTHONPATH=$PYTHONPATH:$main_dir/drio.py
cat $main_dir/downstream/diff_allele_freq/100.vcf |\
  $main_dir/downstream/diff_allele_freq/diff_allele_freq.py \
  $main_dir/downstream/diff_allele_freq/groups.tsv \
  $MIN_STD $MIN_NUM_ALL_PER_GROUP
