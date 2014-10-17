#!/bin/bash
# vim :set ft=2 expandtab:
#
# This script setups a working environment in the current
# directory to run the cnv pipeline from it.
#
# This is how a cnv project should look like:
#
# bam/
# common.sh
# genome/
# step0-rmdups/
# step1-masking/
# step2-read-prep/
# step3-mapping/
# step4-mr-canavar/
# step5-filter-m1/
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Bootstraping structure ..."
cp -r $SRC_DIR/* .
rm -f bootstrap.sh

if [ "$1" != "test" ]; then
(
cat << EOF
# Change this
fasta="`pwd`/genome/XXXXXXXXX"
id="XXXXXXXXXXX"
bam="`pwd`/bam/XXXXXXXXXXXXX"
src="$SRC_DIR/../"

# Do not change
fasta_masked="`pwd`/step1-masking/masked.fa"
fasta_masked_pad="`pwd`/step1-masking/masked.fa.pad"
step1_sh="\$src/step1-masking/run-pipe.sh"
step2_sh="\$src/step2-read-prep/run-pipe.sh"
step3_sh="\$src/step3-mapping/run-pipe.sh"
step4_sh="\$src/step4-mr-canavar/run-pipe.sh"
csv2submit="\$src/csv2submit.py"
EOF
) > common.sh
    echo "Good, edit ./common.sh now."
##############
else
#############
# small.rhemac2.100000.fa
cp $SRC_DIR/../../data/small_rhesus/* ./genome/
curl -L http://cl.ly/013z0Z0q031w/download/24898.300k.bam > bam/24898.300k.bam
curl -L http://cl.ly/173J002G0r2I/download/24898.300k.bam.bai > bam/24898.300k.bam.bai
cp ./genome/*.bed ./step1-masking/chrm_info.bed

(
cat << EOF
fasta="`pwd`/genome/small.rhemac2.100000.fa"
id="foo_test"
bam="`pwd`/bam/24898.300k.bam"
src="$SRC_DIR/../"

fasta_masked="`pwd`/step1-masking/masked.fa"
step0_sh="\$src/step0-rmdups/run-pipe.sh"
step1_sh="\$src/step1-masking/run-pipe.sh"
step2_sh="\$src/step2-read-prep/run-pipe.sh"
step3_sh="\$src/step3-mapping/run-pipe.sh"
step4_sh="\$src/step4-mr-canavar/run-pipe.sh"
csv2submit="\$src/csv2submit.py"

EOF
) > common.sh
fi
