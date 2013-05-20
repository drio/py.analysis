#!/bin/bash

NLINES=10000000
HEAD="head -$NLINES"
#HEAD="cat -"

cat ../locations.bed | $HEAD   > loc.bed
gzip -cd ../30158.bam.depth.txt.gz | $HEAD > depth.txt

cat depth.txt | /stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/add_rdp/add_missing_locations.py loc.bed  > new.bed

ls -lachd depth.txt
ls -lacd new.bed
ls -lacd loc.bed

awk '{print $1"_"$2}' new.bed > n
awk '{print $1"_"$3}' loc.bed > l

diff n l | wc -l
rm -f n l
