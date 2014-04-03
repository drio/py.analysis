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
$ ./run.sh  | sed 's/20000000/50000/g | bash
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
