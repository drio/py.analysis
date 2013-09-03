## Snp calling pipeline (Beta)

NOTE: Please, create a new [issue](https://github.com/drio/py.analysis/issues) if you have
feature requests or find a bug.

### Introduction

This document describes how to use the snp calling pipeline to generate population level snps
for a set of samples.

### Usage (individual snp calls)

Let's assume we have a couple of samples (two bams) we want to get snp calls from. Our final
result should be a population level vcf with the calls from those two samples.

The main entry point for the pipeline is $PIPE/snpcalling-pipe.sh. This script requires certain
arguments in order to run correctly, mainly: path to the directory where the bams are,
the type of dataset (wes or wgs), the path to the reference fasta file and the id of the genome
we want to use for the annotations.

Let's assume we have:

```
$ cd bams
$ ls
35050.bam 35057.bam
```

Now, if we want to generate the jobs to get the snps using rhesus annotations, we would:

```
$ ref="/stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa"
$ bams="`pwd`/bams"
$ DB="MMUL_1.66" # rhesus snpeff db

$ $PIPE/snpcalling-pipe.sh $bams wgs $ref "$DB"
```

When you run the previous command you will get a list of the commands that have to be completed in order to get
the individual snp calls per each sample. If you are happy with the output, you can send them to the cluster
but piping it to your bash:

```
$ $PIPE/snpcalling-pipe.sh $bams wgs $ref "$DB" | bash
```

You can then monitor your jobs in the cluster.

### Usage (Population level vcf)

When the jobs complete, you should have two gzip files with the snp calls in vcf format:

```
$ ls
35057.vcf.anno.gz 35050.vcf.anno.gz
```

Now, to generate the population level vcf we run:

```
$PIPE/merge.snps.sh | bash
```

You should monitor the job. When done, you should get a file named ```merged.vcf.g'```. That's the population
level vcf annotated.

### Downstream analysis ...









