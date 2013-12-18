#!/bin/bash

set -e

FASTA=$1
ONE_TWO=$2
N_THREADS=$3
BAM=$4
SPLIT_NUM=$5
URL=$6
ID=$7

mkdir -p sais
cd sais

bwa aln -t$N_THREADS -b$ONE_TWO $FASTA $BAM > $ONE_TWO.$SPLIT_NUM.sai

signal.py $URL $ID sais

touch done.${SPLIT_NUM}.txt
