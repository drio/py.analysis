### Description

This is a pipeline to compute CNVs on primates.

### Testing with small dataset

Start by bootstraping the *test* project (I assuming you have the pipeline
instaled):

```sh
$ cd /somewhere
$ /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/bootstrap/bootstrap.sh test
Bootstraping structure ...
```

Now it is just a matter of following the steps and executing them in order only when
the previous steps are completed:

```sh
$ cd step1-masking
$ ./run.sh | bash
```

With one minor consideration: When running step2-prep-reads, make sure you change the number of
reads per split since the numbre of reads in the test bam file is very small:

```sh
$ cd step2-read-prep
$ ./run.sh  | sed 's/20000000/50000/g' | bash
```

And continue with the others.

### Working with real data

To work with real data:

```sh
$ mkdir my_project
$ cd my_project
$  ~/dev/py.analysis/pipelines/cnv/src/bootstrap/bootstrap.sh
Bootstraping structure ...
Good, edit ./common.sh now.
```

Now, go ahead and edit `./common.sh` to reflect your data set.

Once done with it, execute the different step as in the testing example.

### A word on disk space usage

This pipeline consumes tons of disk space when used against a mamal genome with 30x coverage samples. 
Here is an example for a 1kg sample:

```sh
$ du -hs *
241G    bam
64K     common.sh
673M    filter
27G     genome
1.3T    step1 # 
11G     step2
9.7T    step3 # <-------- !!!
76G     step4
552M    validations
```

Notice that data in step1 -- particulary alignments --, can be removed once the step completes 
and the masked genome has been generated.

And obiviously, you want to clean up step3 once you have the cnv calls done (filter step). 

In general, you want to clean up everything and keep only the final results and, of course, the 
masked reference genome. That step is particularly time consuming, the good news is you have to 
run it only once per each genome.
