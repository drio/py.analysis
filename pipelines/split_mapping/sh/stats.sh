#!/bin/bash

URL=$1
ID=$2

mkdir -p stats
cd stats

INPUT=`ls ../dups/*.bam`
samtools flagstat $INPUT > stats.txt
signal.py $URL $ID stats
touch done.txt
