#!/bin/bash

TMP_DIR=$1
MEM=$2
INPUT=$3
URL=$4
ID=$5

mkdir -p validate
cd validate

java -Xmx${MEM} -jar $PICARD/ValidateSamFile.jar \
    TMP_DIR=$TMP_DIR \
    INPUT=$INPUT \
    OUTPUT=validate.txt

signal.py $URL $ID validate
touch done.txt



