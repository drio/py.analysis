#!/bin/bash

set -e

BAM=$1
URL=$2
ID=$3

mkdir -p fastqc
cd fastqc
fastqc -o . $BAM

signal.py $URL $ID fastqc

touch done.txt
