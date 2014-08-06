### Description

This document's aim is to describe in detail the status of the efforts towards
building a cnv pipeline. Once the software is completed, it will be validated
using a set of highly analyzed human samples.

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

To test the quality of our calls, plot the copy number calls between one of
the truth set samples and our calls. Notice we have to process both calls since
the window sizes are different. Check ```validation/normalize_windows.py``` for
details on how this is done.

With both call sets in the same windows, we go ahead and dotplot each cnv value.
Notice we collapse adjacent windows to step up the plotting process:

![](http://is04607.com/primates/cnv/dp_calls_vs_truth.png)

Another plot that help us in the validation process is the distribution of
cnv values across the control regions (as defined by mrcanavar):

![](http://is04607.com/primates/cnv/dist_cnv_values_across_control.png)

