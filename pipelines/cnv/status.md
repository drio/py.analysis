### Description

This document's aim is to describe in detail the status of the efforts towards
building a cnv pipeline. Once the software is completed, it will be validated
using a set of highly analyzed human samples.

### Building Blocks

![](https://raw.githubusercontent.com/drio/py.analysis/master/pipelines/cnv/src/schema/cnv.png)

### Step1: Preparing the reference genome. (masked.fa)

In this step we want to mask all the regions of the genome that are over-represented.
To do so, we map the genome against itself by first enumerating adjacent kmers from the
fasta file. We then map those against the genome and identify the locations where we
have been able to align more than one read.

The reference used for validation is: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz

Validation TIP: It would be useful to make sure that two pipeline implementation
mask exactly the same regions. We have done so with two chrms but we could also do
it for the full genome.

### Step2: Prepare sample reads

In this step we iterate over all the reads for the sample (NA12878) and generate two
reads per each 100bp read. The 10 first and last bp of the reads are discarded. We also
split the read dataset to compute alignments in parallel.

### Step3: Map reads

We now map (mrsfast) those reads against the reference genome (masked.fa).

### Step4: Compute genome windows, read depths and call copy number.

We use a modified version of the input genome (masked.fa). We will called masked.fa.pad.

We have done a few things with this new reference:

1. An extra 36bp has been added to the ends of the gaps.
2. The gaps and repeat mask tracks have been used to further mask
those locations in the genome.

```sh
$ curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | \
  gzip -cd - | sed 's/chr//g' | cut -f 6-8 > rMask.bed
$ curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz  | \
  gzip -cd - | sed 's/chr//g' | cut -f 2-4 > trf.bed
```

Finally now use mrcanavar to compute genomic windows, read depths and finally
call copy number on those.

### Validation

Previous attempts on validating the pipeline involved running the pipeline against a 1kb human 
dataset. To improve the validation and since the raw reads for the French sample within the 
truth set are available, we decided to use those.

Once we have results from the pipeline, we compute the overlaps between CNV events called 
from canavar and the events called in the truth set. The script performing this task can 
be found [here](https://github.com/drio/py.analysis/blob/master/pipelines/cnv/src/validation/v2-run-intersect.sh).

These are the results:

```
truth/originals/Homo_sapiens-French_HGDP00521_1000bp_simple_per_1kb.HM.bedgraph.bed (-1, 3)
pipeline.calls.bed (-1, 3, 100)

overlap bp_overlap num_events_a bp_events_a num_events_b bp_events_b num_overlaps_a_b bp_overalps_a_b
--
A (Truth) ∩ B (Pipeline calls): .63 .80 38999 217112153 191698 1459170215 24720 174655789
B (Pipeline calls) ∩ A (Truth): .12 .26 191698 1459170215 38999 217112153 24339 387943058
```

It is important to notice that the posible cnv values for the truth set calls are [0 .. 10].
In our calls, the range of possible values for the cnv values are: [0.. 100000]. High values for 
the cnv calls may be artifact. The script only accounts events that are smaller than a certain
value (*<100*).

Notice the interval values at the end of the first two lines. Those are the min and max values 
that we use as threshold to determine what calls we include when computing the intersections.
-1 means we are ignoring deletions and only accounting for duplications.

The script computes the overlap between A/B and B/A. A being the cnv events from the truth set
and B the events called by our pipeline.

The results show a 63% overlap between A and B (80 if we look at bp level). When we compute the overlaps
the other way around we get a 12% overlap (.26% if we look at bp level).



