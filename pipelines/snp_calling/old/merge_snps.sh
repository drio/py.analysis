#!/bin/bash
#
# merge all snps and annotate the merged file
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

for f in `find . -name "*.vcf*.gz"`
do
  input="$input $f"
done

out=merged.vcf.gz
cmd="vcf-merge $input | bgzip -c > $out && tabix -p vcf $out"
echo "$cmd" | submit -s "merge.anno.ex" -m15G -c4
