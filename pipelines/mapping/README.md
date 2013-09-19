# Basic mapping pipeline

## When to use this

Use this pipeline to compute aliments for Illumina PE/MP data.

The steps are:

  1. sais
  2. sampe
  3. sam->bam + sorting
  4. coordinate sorting
  5. remove temporary files

The entry point is ```jobs_for_mapping.sh```. That script will generate the jobs files
including the necessary dependencies to compute the whole mapping process

## Quickstart

Here you have an example when running the pipeline in a single machine.

```sh
  $ export PATH=$PATH:/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/mapping
  $ export PATH=$PATH:/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines
  $ jobs_for_mapping.sh -b /Users/drio/dev/bam.examples/phix.bam \
    -f /Users/drio/dev/genomes/phix.fa \
    -o FOOO -r 1G -m /tmp -s 11111 -t 4 | awk -F"\t" '{print $5}'  | grep -v cmd
```

```jobs_for_mapping.sh``` generates a table with all the cmds and its dependencies necessary to compute the
mappings.

Run jobs_for_mapping.sh without parameters to get an explanation of what each flag does.

Before going any further you should run this example to make sure your environment is properly set up. You
are probably going to do it in ardmore. Here you have some important locations:

If working in ardmore, Phix bam and reference can be found here:

```sh
/stornext/snfs6/rogers/drio_scratch/playground/bam.examples/phix.bam
```

```sh
/stornext/snfs6/rogers/drio_scratch/genomes/phix.fa
```

So, you should run something like this:

```sh
$ export PATH=$PATH:/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/mapping
$ export PATH=$PATH:/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines
$ cd somewhere/in/your/directories
$ jobs_for_mapping.sh -b /stornext/snfs6/rogers/drio_scratch/playground/bam.examples/phix.bam \
  -f /stornext/snfs6/rogers/drio_scratch/genomes/phix.fa -o FOOO -r 1G -m /tmp -s 1111 -t 2 |\
  awk -F"\t" '{print $5}'  | grep -v cmd | bash
....
$ ls -l *.bam
-rw-r--r-- 1 deiros rogers 3521971 Sep 19 09:07 FOOO.sorted.dups.bam
```

## Other features

### Computing certain parts of the pipeline

If you have computed part of the jobs, you can drop some of them with some unix magic.
Let's say you don't want to recompute the sais:

```sh
  $ jobs_for_mapping.sh -b /stornext/snfs6/rogers/drio_scratch/playground/bam.examples/phix.bam \
    -f /stornext/snfs6/rogers/drio_scratch/genomes/phix.fa -o FOOO -r 1G -m /tmp -s 1111 -t 2 |\
    awk -F"\t" '{print $5}'  | grep -v cmd |\
    grep -vP ".sai$" | ruby -ne 'puts $_.sub(/\d+,\d+/, "-")' | bash
```

The regex part of the pipeline is necessary because the id of the jobs are
pseudo-randomly generated.

### Running in a cluster

If you are planning to run this in a cluster (pbs), use ```cmds2submit.py``` to generate
the actual submit commands:

```sh
  $ jobs_for_mapping.sh -b /stornext/snfs6/rogers/drio_scratch/playground/bam.examples/phix.bam \
    -f /stornext/snfs6/rogers/drio_scratch/genomes/phix.fa -o FOOO -r 1G -m /tmp -s 1111 -t 2 |\
    awk -F"\t" '{print $5}'  | grep -v cmd | cmds2submit.py -
```

At this point is probably a good idea to run the canonical example (phix) in cluster mode so you 
are sure things are working properly.

