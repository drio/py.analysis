This tool adds real coverage information per each sample to a merged vcf file.

We have to generate first the bed file with the coverage per each sample:

chrm_coor cov_sample1 cov_sample2 ... cov_sampleN
Example: chr1_123123123 12 34 ... 32

You can use the join_coverate.sh as a template. This will generate the
locations.bed and the joined coverage file

Once that is done we can run this tool:

$ cat coverage.txt | add_rdp.py merged.vcf.gz | gzip -c - > merged.rdp.gz


Notice these:

num(merged.vcf) == num_lines(locations.bed) > num_lines(join.bed)

The reason for the number of lines in join.bed being less than the number of snps
is because when looking for depth coverage, we require the reads to have at least
MAPQ = 1. If all the reads in a locus have 0 MAPQ, samtools depth will not report
that locus.
