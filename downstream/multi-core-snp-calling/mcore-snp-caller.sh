#!/bin/bash
#
set -e

snp_eff_dir=/stornext/snfs6/rogers/drio_scratch/local/snpeff/snpEff_3_0a
snp_eff_db="MMUL_1.66"


error() {
  echo "$1"
  cat <<-EOF
Usage:
  $ ./experimental-caller.sh <bam file> <path to ref> [num of concurrent processes]
  If num of conc. is not provided we will use the number of cores available.

Example:
  $ ./experimental-caller.sh 9498-01A.bam /stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa 10

 NOTES:
  1. You probably want to get a big machine: msub -I -q gac -l nodes=1:ppn=16,mem=50Gb -l feature=bigmem -d \`pwd\`
  2. Indivials instances don't use too much memory, but if you run 32 instances you can eat up quite a lot of memory.
     Keep that in mind when requesting memory for your interactive job. Each instance will eat up about 6G for a human sample.
EOF
  exit 1
}

bam=$1
ref_fasta=$2
n_concurrent=$3

if [ ! -f "$bam" ] || [ ! -f $bam.bai ]
then
  error "need indexed bam"
fi

if [ ".$ref_fasta" == "." ] || [ ! -f $ref_fasta ]
then
  error "need ref fasta"
fi

if [ ".$n_concurrent" == "." ]
then
  n_concurrent=`cat /proc/cpuinfo  | grep proc | wc -l`
fi

# Compute likelihoods per each possible genotype
# http://massgenomics.org/2012/03/5-things-to-know-about-samtools-mpileup.html
# -S   output per-sample strand bias P-value in BCF (require -g/-u)
# -u   uncompress BCF output
# -f   fasta file (indexed)
# -q1  skip alignments with mapQ smaller than 1
id=`basename $bam | sed 's/.bam//'`

# Compute concurrently the snps per each chrm
rm -f $id.*.anno.vcf
tmp_chrm_vcfs=$id.{}.anno.vcf
samtools view -H $bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | \
  xargs -t -I {} -n 1 -P $n_concurrent sh -c \
  "samtools mpileup -q1 -S -uf $ref_fasta -r {} $bam | \
    bcftools view -vcg - | \
    java -Xmx4G -jar $snp_eff_dir/snpEff.jar eff -c $snp_eff_dir/snpEff.config -v $snp_eff_db -noStats > $tmp_chrm_vcfs"

# merge the results
tmp_header=$id.header.tmp
tmp_snps=$id.tmp
rm -f $tmp_header $tmp_snps
find . -name "$id.*.anno.vcf" | head -1 | xargs -I {} cat {} | grep "^#" > $tmp_header
find . -name "$id.*.anno.vcf" | xargs -I {} cat {} | grep -v "^#" >> $tmp_snps

# Generate final results (vcf)
final_output=$id.anno.vcf.gz
cat $tmp_header $tmp_snps | bgzip -c > $final_output
tabix -p vcf $final_output

# Clean up
rm -f $tmp_header $tmp_snps $id.*.anno.vcf.gz
