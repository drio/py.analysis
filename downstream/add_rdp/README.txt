This tool adds real coverage information per each sample to a merged vcf file.

We have to generate first the bed file with the coverage per each sample:

chrm_coor cov_sample1 cov_sample2 ... cov_sampleN
Example: chr1_123123123 12 34 ... 32

You can use the all.sh as a template. This will:

1. Generate the join coverage for all the samples
2. run the add_rdp tool to add the rdp to the vcf file
3. optional: filter the final answer based on snp quality and the rdp.

Notice these:

num(merged.vcf) == num_lines(locations.bed) > num_lines(join.bed)

The reason for the number of lines in join.bed being less than the number of snps
is because when looking for depth coverage, we require the reads to have at least
MAPQ = X, where X is > 0. If all the reads in a locus have a MAPQ < X, samtools depth
will not report that locus.
