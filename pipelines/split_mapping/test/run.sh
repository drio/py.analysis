#!/bin/bash

[ $(which sapi.py) ] || (echo "sapi.py not found"; exit 1)

curl -L http://cl.ly/4638453t2Q0E/phix_dataset.tar.bz2 |  tar -jx


bam="`pwd`/phix/phix.bam"
fa="`pwd`/phix/phix.fa"
sapi.py -b $bam fastqc | bash
sapi.py -b $bam  -n 10000 splits | bash
sapi.py -b $bam -f $fa sais | bash
sapi.py -b $bam  -f $fa sampe | bash
sapi.py merge | bash
sapi.py dups | bash
sapi.py stats | bash
