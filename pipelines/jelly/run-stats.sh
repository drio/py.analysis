#!/bin/bash

echo "perl ./assemblathon_stats.pl  ../reference/atlas.fasta > input.stats.txt" | submit -s input
echo "perl ./assemblathon_stats.pl  ../jelly.out.fasta > jelly.stats.txt" | submit -s jelly
