PIPE_DIR=/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv
cp small_rhesus/small.rhemac2.100000.sizes.bed ./chrm_info.bed

mrsfast --index small_rhesus/small.rhemac2.100000.fa --ws 12 small_rhesus/small.rhemac2.100000.fa.index

$PIPE_DIR/src/step1-masking/makeIntervalsBED_v2.pl chr1 100000 36 5 > chr1_k36_step5_intervals.bed

fastaFromBed -fi small_rhesus/small.rhemac2.100000.fa -bed chr1_k36_step5_intervals.bed -fo chr1.fa

mrsfast --search small_rhesus/small.rhemac2.100000.fa --seq chr1.fa -o /dev/stdout -e 2 | cut -f1 | gzip -c > chr1_k36_step5_intervals.map.gz

gzip -cd chr1_k36_step5_intervals.map.gz | /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/step1-masking/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl chr1.fa chr1 > counts.chr1.bed

cat counts* | grep -v \# | grep -v -P \^$ >> counts.bed

maskFastaFromBed -fi small_rhesus/small.rhemac2.100000.fa -bed counts.bed -fo masked.fa
cat counts.bed | R CMD BATCH /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/step1-masking/makeHistogramCounts.R

