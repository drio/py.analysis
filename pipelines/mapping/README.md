# Basic mapping pipeline

## What
Use this pipeline to compute aliments for Illumina PE/MP data.

The steps are:

  1. sais
  2. sampe
  3. sam->bam + sorting
  4. coordinate sorting
  5. remove temporary files

The entry point is ```jobs_for_mapping.sh```. That script will generate the jobs files
including the necessary dependencies to compute the whole mapping process

## Basic Example

Here you have an example when running the pipeline in a single machine.

```sh
  $ export PATH=$PATH:/path/to/jobs_for_mapping_dir # In ardmore: /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/mapping
  $ jobs_for_mapping.sh -b /Users/drio/dev/bam.examples/phix.bam \
    -f /Users/drio/dev/genomes/phix.fa \
    -o FOOO -r 1G -m /tmp -s 11111 -t 4 | awk -F\t '{print $5}'  | grep -v cmd
```

```jobs_for_mapping.sh``` generates a table with all the cmds and its dependencies necessary to compute the
mappings.

Run jobs_for_mapping.sh without parameters to get an explanation of what each flag does.

## Other features

### Computing certain parts of the pipeline

If you have computed part of the jobs, you can drop some of them with some unix magic.
Let's say you don't want to recompute the sais:

```sh
  $ jobs_for_mapping.sh /Users/drio/dev/bam.examples/phix.bam \
      /Users/drio/dev/genomes/phix.fa FOOOOO 1G /tmp 4 | \
      grep -vP ".sai$" | ruby -ne 'puts $_.sub(/\d+,\d+/, "-")' | cmds2submit.py -
```

The regex part of the pipeline is necessary because the id of the jobs are
pseudo-randomly generated.

### Running in a cluster

If you are planning to run this in a cluster (pbs), use ```cmds2submit.py``` to generate
the actual submit commands:

```sh
  $ jobs_for_mapping.sh /Users/drio/dev/bam.examples/phix.bam /Users/drio/dev/genomes/phix.fa FOOOOO 1G /tmp 4  | ./cmds2submit.py -
```




