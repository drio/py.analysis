Intro
=====

This is a simple pipeline to use mapping stats to determine sample swaps.
What we do is to extract a random subset of reads from the input fastq, 
map those against are reference and generate stats on the results.


Usage
-----

First, create a directory where you'll be working on, and within it, create
a `src` directory that points to the source of the pipeline

```sh
$ mkdir my_project 
$ cd my_project
$ ln -s /stornext/snfs7/rogers/drio_scratch/dev/py.analysis/downstream/sample-swap-check ./src
```

Now, let's say you have some fastq file you want to process with the pipeline and 
also a reference file you want to use for mapping (NOTE: this reference has to 
be processed with bwa to generate the indexes).

```sh
$ ls -aclhd data/ok_baboon/s_2_1_sequence.txt.bz2
lrwxrwxrwx 1 deiros rogers 226 Nov 11 14:24 data/ok_baboon/s_2_1_sequence.txt.bz2 -> /stornext/snfs7/rogers/baboondiversity_fastq/./restores/20140619/stornext/snfs0/next-gen/Illumina/Instruments/700142/101122_SN142_0253_A805CEABXX/Data/Intensities/BaseCalls/GERALD_07-12-2010_p-illumina.4/s_2_1_sequence.txt.bz2
$ ls -aclhd genomes/rhemac2.indian_macaque_no_phix.fixed.fa
-rw-r--r-- 1 deiros rogers 3.0G Jun 14 16:24 genomes/rhemac2.indian_macaque_no_phix.fixed.fa
```

Great. Now we can test the pipeline with a small set of reads (2500):

```sh
$ bzip2 -cd data/ok_baboon/s_2_1_sequence.txt.bz2 | head -10000 | ./src/run-pipe.sh genomes/rhemac2.indian_macaque_no_phix.fixed.fa foo 1000
[bwa_aln] 17bp reads: max_diff = 2
[bwa_aln] 38bp reads: max_diff = 3
[bwa_aln] 64bp reads: max_diff = 4
[bwa_aln] 93bp reads: max_diff = 5
[bwa_aln] 124bp reads: max_diff = 6
[bwa_aln] 157bp reads: max_diff = 7
[bwa_aln] 190bp reads: max_diff = 8
[bwa_aln] 225bp reads: max_diff = 9
[bwa_aln_core] calculate SA coordinate... 7.00 sec
[bwa_aln_core] write to the disk... 0.00 sec
[bwa_aln_core] 1000 sequences have been processed.
[main] Version: 0.6.2-r126
[main] CMD: bwa aln -t4 genomes/rhemac2.indian_macaque_no_phix.fixed.fa foo.aORfKp2XbJs=.fq
[main] Real time: 48.290 sec; CPU: 25.102 sec
[bwa_aln_core] convert to sequence coordinate... 32.84 sec
[bwa_aln_core] refine gapped alignments... 5.53 sec
[bwa_aln_core] print alignments... [samopen] SAM header is present: 44 sequences.
0.02 sec
[bwa_aln_core] 1000 sequences have been processed.
[main] Version: 0.6.2-r126
[main] CMD: bwa samse genomes/rhemac2.indian_macaque_no_phix.fixed.fa /dev/fd/63 foo.aORfKp2XbJs=.fq
[main] Real time: 150.396 sec; CPU: 38.766 sec
```

The only thing we are doing here is to extract just 10000 lines from our
fastq file and feeding that into the pipeline for processing. This takes
about 2 minutes in my cluster. The arguments for the pipeline are these:

```
$ ./src/run-pipe.sh
ERROR:  Need sample name as second parameter
Usage: run-pipe.sh <ref.fa> <sample_name> <n_reads_to_sample>
```

Now when you are done you should see a file called `stats.foo.1000.XXXXXXX.txt` 
where XXXXXXX is a weird unique string, foo is the sample name you passed and 
1000 is the number of reads you wanted to extract from the original fastq data.

The contents of this file look like this:

```
$ cat stats.foo.1000.aORfKp2XbJs\=.txt
1000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
456 + 0 mapped (45.60%:-nan%)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (-nan%:-nan%)
0 + 0 with itself and mate mapped
0 + 0 singletons (-nan%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Using the pipeline do detect sample swaps
-----------------------------------------

If you want to detect sample swaps, what you can do is to run the
"conflictive" dataset (fastq) against the expected reference and
also against other genomes (typically human). You probably want 
to run it multiple times with a reasonable amount of reads (perhaps
10 million). This will mean multiple stats outputs.

Then you can compare the stats output. If you are expecting, let's
say, baboon, you should have a consistently higher number of mapped
reads to the baboon reference and less reads mapping to the human
reference. Or the other way around if there is a swap.




