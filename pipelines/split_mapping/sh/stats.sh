#!/bin/bash

mkdir -p stats
cd stats

INPUT=`ls ../dups/*.bam`
samtools flagstat $INPUT > stats.txt
