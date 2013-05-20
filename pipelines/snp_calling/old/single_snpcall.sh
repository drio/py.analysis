#!/bin/bash
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SRC_DIR/../common.sh

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
#
# NOTE: We are not doing any read depth filtering. A couple of things:
# 1. A very few numbers of snps will be generated with low coverage
# 2. The number of snps with low read depth will have a low snp quality
# So, you can take care of those very easily downstream by looking at the
# snp quality.
#
samtools mpileup -q$MIN_MAP_QUALITY -S -uf $ref_fasta $ID.bam | bcftools view -bvcg - > $ID.bcf

# Do the actual calling + annotate the snps (if DB provided)
if [ ".$snp_eff_db" == "." ];then
  bcftools view $ID.bcf | bgzip -c > $ID.vcf.gz
  out="$ID.vcf.gz"
else
  bcftools view $ID.bcf | \
    java -Xmx4G -jar $SNP_EFF_DIR/snpEff.jar eff -c $SNP_EFF_DIR/snpEff.config -v $snp_eff_db -noStats |\
    bgzip -c > $ID.vcf.anno.gz
  out="$ID.vcf.anno.gz"
fi

in=\$out
# Index the file
tabix -p vcf \$in

# Clean up
rm -f $ID.bcf
EOF
