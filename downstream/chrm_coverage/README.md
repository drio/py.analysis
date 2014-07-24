### Intro

This is a simple pipeline to generate per read coverage plots on
bam files.

Given a bam file with alignments, it will compute the per base coverage
`coverage.bed.gz` and generate per chrm coverage plots.

### Usage

```sh
$ mkdir dir
$ cd dir
$ ln -s /path/to/bam/to/analyze ./input.bam
$ ln -s /stornext/snfs7/rogers/drio_scratch/dev/py.analysis/downstream/chrm_coverage/Makefile
$ make plot.Chr1.png
```

The name of the chrm (`Chr1` in this case), has to match the name of used when
generating the alignments. Check the header of the bam if you are not sure.

In addition, if you are planning to compute the plots for all the chromosomes
in parallel, make sure you compute the coverage first by running: `make
coverage.bed.gz`.

By default, the genomic window size used for plotting is 10kb. A downsampling
is necessary since the input space (bp region) is much bigger than the plotting
space (pixels).
