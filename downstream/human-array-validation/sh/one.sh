#!/bin/bash

#/stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/human-array-validation/test_data/probes.fq \
#/stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/human-array-validation/test_data/18277.bam \
#/stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/human-array-validation/test_data/affy6.fixed.txt \

./src/compute_a_sample.sh \
./probes/affy6.txt.fixed.readified.fq \
bams/18277.bam \
18277 \
/stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa \
./probes/affy6.txt.fixed \
affy6
