
rand_string() {
  [ `which md5 2>/dev/null` ] && MD5=md5
  [ `which gm5sum 2>/dev/null` ] && MD5=gm5sum
  echo "$RANDOM" | $MD5 | awk '{print $1}'
}

#MIN_COV=3
#MAX_COV=1000
# valid databases: MMUL_1.66, GRCh37.66
SNP_EFF_DIR=/stornext/snfs6/rogers/drio_scratch/local/snpeff/snpEff_3_0a
MIN_MAP_QUALITY=10
#MIN_COV=3
#MAX_COV=1000
#SNP_EFF_DB="MMUL_1.66"
#REF_FASTA=/stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa



