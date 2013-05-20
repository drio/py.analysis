#!/bin/bash

GS="2900000000"; for i in *.gz;do echo -ne "$i: "; ns=`zcat $i | grep -v "#"| grep -v INDEL | wc -l`; echo "scale=10; ($ns/$GS)*1000"| bc ;done  > $HOME/snps.exon.txt
