#!/bin/bash
#
# merge all snps and annotate the merged file
#

for f in `find . -name "*.vcf*.gz" | grep -v merged`
do
  input="$input $f"
done

snp_eff_dir=/stornext/snfs6/rogers/drio_scratch/local/snpeff/latest
#out=merged.anno.rdp.vcf.gz
out=merged.vcf.gz
#anno="java -Xmx4G -jar $snp_eff_dir/snpEff.jar eff -c $snp_eff_dir/snpEff.config -v MMUL_1.66 -noStats"
cmd="vcf-merge $input | bgzip -c > $out && tabix -p vcf $out"
echo "$cmd" | submit -s "merge.anno.ex" -m15G -c4
