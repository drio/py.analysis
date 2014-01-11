### Objective

We want to determine how many genes would be captured effectively when using
rhesus DNA on a human chip.

### Pipeline steps (first approach)

1. Download list of annotated genes from ensembl (gtf)
2. Find read depth for those coordinates
3. Generate descriptive statistics on the exon locations for a particular sample
(bam; reads mapped to rhesus).
4. Determine what exons are filtered out and which ones pass the filtering.

### Filtering criteria

Very simple: find the average read depth for all the
exons.  Then use that average to create an interval of valid read depth.  Any
exon with a read depth within the interval will pass the filter.

The interval will be generated like this:

```
average = compute average RD on all exons for that sample
interval_size = average - 20
min = average - interval_size
max = average + interval_size
interval = (min, max)
```

### Using the pipeline (single sample)

```sh
$ MAKE=/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/capture_efficiency/Makefile
$ make -f $MAKE ID=18277 BAM=bams/18277.bam
```

### Using the pipeline (multiple samples)

```sh
# Notice here I am using submit which takes care of sending jobs to your cluster
$ (for b in bams/*.bam; do id=`basename $b | sed 's/.bam//'`; submit -s ${id}_crv "make -f $MAKE ID=$id BAM=$b"; echo "#"; done) | bash
```

### Others

You can use ```src/means.py``` to get the mean of ane exon for all the samples. This
may be useful to perform further filtering.
