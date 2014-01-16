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

### Population level filtering

Once we have filtered out the targets (exons) for the individuals samples we can apply
a second level of filtering. From the same directory where we generated the individual
results:

```sh
$ make -f /path/to/Makefile all.out.10.bed all.pass.10.bed
```

Notice the software will use ```.10.``` in the filtering. In this case, we have to have
at least 10 samples for which that exon passed the first level filtering.

```sh
$ head all.out.10.bed
9       68868671        68868838        RUFY2   9       RUFY2-203       6
5       157418332       157418453       MMU.2733        4       MMU.2733-201    8
4       105767394       105767552       ARID1B  9       ARID1B-203      4
20      2904690 2904716 PRSS21  2       PRSS21-201      9
2       45808430        45808574        MMU.10529       4       MMU.10529-201   9
...
$ head all.pass.10.bed
13      9950891 9950908 TAF1B   1       TAF1B-201       47
5       106107700       106107798       MMU.4610        12      MMU.4610-205    38
3       176871164       176871258       MMU.13875       2       MMU.13875-201   49
3       176871164       176871258       MMU.13875       2       MMU.13875-203   49
3       176871164       176871258       MMU.13875       2       MMU.13875-202   49
```
