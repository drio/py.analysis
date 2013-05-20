#!/bin/bash

ORIGIN="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/../../py.analysis
DST="/stornext/snfs6/rogers/drio_scratch/"

rsync -avz \
--exclude=".git"  --exclude="data/merged.wes.500k.vcf.gz" \
--exclude="py.analysis/downstream/human-array-validation/test_data/18277.bam" \
--progress $ORIGIN ardmore:$DST
