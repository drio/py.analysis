#!/bin/bash

TMP_DIR=$1
MEM=$2
INPUT=$3

mkdir -p validate
cd validate

java -Xmx${MEM}g -jar $PICARD/ValidateSamFile.jar \
    TMP_DIR=$TMP_DIR \
    INPUT=$INPUT \
    OUTPUT=validate.txt




