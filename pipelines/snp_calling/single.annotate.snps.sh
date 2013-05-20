#!/bin/bash
#
# annotate and index single vcfs
#
snp_eff_dir=/stornext/snfs6/rogers/drio_scratch/local/snpeff/snpEff_3_0a

for f in `find bak -name "*snps.all.gz"`
do
  id=`echo $f | ruby -ne 'puts $_.match(/(\d+)/)[1]'`
  out=$id.anno.vcf.gz
  cmd="gzip -cd $f |\
    java -Xmx4G -jar $snp_eff_dir/snpEff.jar eff -c $snp_eff_dir/snpEff.config -v MMUL_1.66 -noStats |\
    bgzip -c > $out && tabix -p vcf $out"
  echo "$cmd" | submit -s "an.ex.$id"
done
