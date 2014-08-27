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
pipeline.calls.bed (-1, 3)

overlap bp_overlap num_events_a bp_events_a num_events_b bp_events_b num_overlaps_a_b bp_overalps_a_b
--
A ∩ B: .65 .90 38999 217112153 191755 1685806730 25706 195511042
B ∩ A: .12 .34 191755 1685806730 38999 217112153 24515 5
```

It is important to notice that the posible cnv values for the truth set calls are [0 .. 10].
In our calls, the range of possible values for the cnv values are: [0.. 100000]. High values for 
the cnv calls may be artifact. The script only accounts events that are smaller than a certain
value (*<100*).

The script computes the overlap between A and B and B and A. A being the cnv events from the truth set
and B our pipeline calls.

The results show a 65% overlap between A and B (90 if we look at bp). When we compute the overlaps
the other way around we get a 12% overlap (.34% if we look at bp level).



