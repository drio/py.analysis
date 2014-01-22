#!/bin/bash
#
set -e

chrm="chr19"
size=1000000
mrsfast=mrsfast/bin/mrsfast
fasta_genome=data/chr19.1M.fa

perl src/step1/makeIntervalsBED_v2.pl $chrm $size > intervals.bed
fastaFromBed -fi $fasta_genome -bed intervals.bed -fo kmers.fa
$mrsfast --index $fasta_genome --ws 12 ${fasta_genome}.index
$mrsfast --search $fasta_genome --seq kmers.fa -o out.map -e 2

src/step1/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl kmers.fa out.map $chrm > counts.bed
#cat $fasta_genome | src/step1/calculateAssemblyStats_v2.pl > stats.txt
cat counts.bed | R CMD BATCH src/step1/makeHistogramCounts.R
maskFastaFromBed -fi $fasta_genome -bed counts.bed -fo masked.fa

echo "rm -rf out.map* counts.txt stats.txt  *.Rout *.png *.pdf *.bed *.txt *.fa"
