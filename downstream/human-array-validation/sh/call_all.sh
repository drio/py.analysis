#!/bin/bash

# Compute all the calls per each sample
i=0
MIN_MAPQ=20
PLUS_MINUS_HITS=5
for f in *counts*.txt
do
  i=$[$i+1]
  p_seed=`echo $f | ruby -ne 'puts $_.match(/counts\.([\w\.-]+)\.fixed\.txt/)[1]'`
  sam=$p_seed.fixed.sam

  echo -ne "$f "
  if [ ! -f $f.calls ];then
    ./src/caller.py $sam $f $MIN_MAPQ $PLUS_MINUS_HITS > $f.calls
    if [ $i -eq 14 ];then
      wait
      i=0
    fi
  fi
done
wait

# Find, per each chip (3), what probes pass for all the samples
for s in Omni affy6 exon
do
  echo -ne "$s: "
  cat *$s*.calls |  grep pass | awk '{print $1}' | sort | uniq -c | grep 55 | wc -l
done
