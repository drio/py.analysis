## Split mapping pipeline

### What is this?

The amount of reads generated in a single illumina lane is increasing constantly. Mapping all those reads is quite straight forward thanks to aligners like BWA,
doing so in an efficient way (both time and space) is more complicated. This pipeline tries to help with that.


### Quickstart

NOTE: For hgsc people (Trivial to change for general use)

#### (1) Setup your environment:

```$ source /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/load_env```

#### (2) Run the pipeline now:

```
$ sapi.py
usage: sapi.py [-h] [-n NUM_READS] [-b BAM] [-f FASTA] [-t N_THREADS]
               [--scheduler {single,pbs}]
               step
sapi.py: error: too few arguments
```

#### (3) Run the pipeline against the test dataset (phix):

```
$ mkdir test; cd test
$ /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/split_mapping/test/run_phix.sh
```

You are going to see a lot of output and after a minute or so you should have an output like this:

```
$ ls
dups  fastqc  merge  phix  sais  sampe  splits  stats
```

Your final bam some stats will be found here:

```
$ ls -aclhd dups/merged.sorted.dups.bam
-rw-r--r-- 1 deiros rogers 4.3M Nov 19 11:46 dups/merged.sorted.dups.bam

$ cat stats/stats.txt
100000 + 0 in total (QC-passed reads + QC-failed reads)
41125 + 0 duplicates
62587 + 0 mapped (62.59%:-nan%)
100000 + 0 paired in sequencing
50000 + 0 read1
50000 + 0 read2
14732 + 0 properly paired (14.73%:-nan%)
40720 + 0 with itself and mate mapped
21867 + 0 singletons (21.87%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

If you got to this point, congrats you are ready to use the pipeline with your data.

### Real datasets

```sh
$ sapi-bootstrap.sh
// edit your config.sh
$ ./run.sh
```

