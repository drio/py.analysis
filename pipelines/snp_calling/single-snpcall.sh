#!/bin/bash
#
PRJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BAMS_DIR=${PRJ_DIR}/bams
#MIN_COV=3
#MAX_COV=1000
SNP_EFF_DIR=/stornext/snfs6/rogers/drio_scratch/local/snpeff/latest
REAL_DP_CLASS="/stornext/snfs6/rogers/drio_scratch/drd.bio.toolbox/java/realDP-vcf"

SNP_EFF_DB="MMUL_1.69"
REF_FASTA=/stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa

help() {
  local msg=$1
  [ ".$msg" != "." ] && echo "UPS!: $msg"
  echo "Usage: single.sh <sample_id> <path_to_bam> <wes_or_wgs> <ref_fasta> <snp_eff_db_id>"
  exit 1
}

ID=$1
bam=$2
wes_or_wgs=$3
ref_fasta=$4
snp_eff_db=$5
[ ".$ID" == "." ] && help
[ ".$bam" == "." ] && help
[ ".$ref_fasta" == "." ] && help
[ ! -f "$ref_fasta" ] && help "I couldn't open the fasta file."
if [ ".$snp_eff_db" == "." ];then
  echo "warning: NO snp effect db selected. We won't annotate!" >&2
fi

cat <<EOF
#!/bin/bash
set -e
# BAM = $bam
ln -s $bam ./$ID.bam
samtools index ./$ID.bam

# Compute likelihoods per each possible genotype
# http://massgenomics.org/2012/03/5-things-to-know-about-samtools-mpileup.html
# -S   output per-sample strand bias P-value in BCF (require -g/-u)
# -u   uncompress BCF output
# -f   fasta file (indexed)
# -q10 skip alignments with mapQ smaller than 10
samtools mpileup -q10 -S -uf $ref_fasta $ID.bam | bcftools view -bvcg - > $ID.bcf

# Do the actual calling + annotate the snps
if [ ".$snp_eff_db" != "." ];then
  bcftools view $ID.bcf | \
    java -Xmx4G -jar $SNP_EFF_DIR/snpEff.jar eff -c $SNP_EFF_DIR/snpEff.config -v $snp_eff_db -noStats |\
    bgzip -c > $ID.vcf.anno.gz
  out="$ID.vcf.anno.gz"
else
  bcftools view $ID.bcf | bgzip -c > $ID.vcf.gz
  out="$ID.vcf.gz"
fi

in=\$out
# Index the file
tabix -p vcf \$in

# RDP
# valid databases: MMUL_1.66, GRCh37.66
#gzip -cd $ID.vcf.anno.gz |\
#  java -Xmx4g -classpath "$PICARD/sam-1.79.jar:$PICARD/picard-1.79.jar:.:$REAL_DP_CLASS" vcfAddCoverage |\
#  bgzip -c > $ID.vcf.anno.rdp.gz

# Filter snps
# bgzip -cd $ID.vcf.anno.gz | filter-vcfs.py "$wes_or_wgs" "sample" | bgzip -c > $ID.vcf.anno.filtered.gz

# Clean up
rm -f $ID.bcf
EOF
