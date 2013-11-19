#!/bin/bash

curl -L http://cl.ly/4638453t2Q0E/phix_dataset.tar.bz2 |  tar -jx

export PATH="$PATH:../"
sapi.py -b $HOME/dev/bam.examples/phix.bam fastqc | bash
sapi.py -b $HOME/dev/bam.examples/phix.bam  -n 10000 splits | bash
sapi.py -b $HOME/dev/bam.examples/phix.bam  -f $HOME/dev/genomes/phix.fa sais | bash
sapi.py -b $HOME/dev/bam.examples/phix.bam  -f $HOME/dev/genomes/phix.fa sampe | bash
sapi.py merge | bash
sapi.py dups | bash
sapi.py stats | bash
