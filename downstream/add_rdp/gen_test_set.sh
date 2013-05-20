#!/bin/bash
#

gzip -cd ../merged.vcf.gz | head -10000 | gzip -c > merged.vcf.gz &
for b in ../*.bam
do
  bn=`basename $b`
  (samtools view -h $b | head -300000 | samtools view -Sb - > $bn; samtools index $bn) &
done
wait
