#!/bin/bash
# Enumerate all the exons per each sample associating the filtering flag
#
for f in *.pass.gz *.out.gz;do
  id=`echo $f | ruby -ane 'puts $_.match(/(\d+)/)[1]'`
  gzip -cd $f | awk -v id=$id '{print id"\t"$1"\t"$2"\t"$3"\t"$15}' | sed 's/Chr//g'
done
