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

### Group level filtering

Here we want to make sure there are enough samples that confirm the hypothesis
that the exon has been captured successfully:

```sh
$ make -f /path/to/Makefile all.out.10.bed all.pass.10.bed
```

That will generate a single all.bed file with all the exons and it will print
in the last column the number of samples that confirm the quality of the exon.
It will also generate one file for the exons that pass and do not pass the
group filtering (In this example I am requiring 10 samples).

```sh
$ head all.bed
chrm    start   end     gene    exon_number     transcript_number       num_samples_pass
13      9950891 9950908 TAF1B   1       TAF1B-201       47
5       106107700       106107798       MMU.4610        12      MMU.4610-205    38
3       176871164       176871258       MMU.13875       2       MMU.13875-201   49
19      21849838        21849902        ZNF492  8       ZNF492-204      0
3       176871164       176871258       MMU.13875       2       MMU.13875-203   49
3       176871164       176871258       MMU.13875       2       MMU.13875-202   49
15      109580808       109580822       MMU.17419       2       MMU.17419-202   45
17      14931782        14931937        NBEA    26      NBEA-201        45
15      92045796        92045857        FRMD3   4       FRMD3-204       43
...

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
...
```
