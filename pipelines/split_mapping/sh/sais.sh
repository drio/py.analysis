#!/bin/bash

FASTA=$1
ONE_TWO=$2
N_THREADS=$3
BAM=$4
SPLIT_NUM=$5

mkdir -p sais
cd sais

bwa aln -t$N_THREADS -b$ONE_TWO $FASTA $BAM > $ONE_TWO.$SPLIT_NUM.sai