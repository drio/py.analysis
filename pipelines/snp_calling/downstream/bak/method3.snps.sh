#!/bin/bash

test_ds="18277.snps.all.gz  31506.snps.all.gz"
sort_mem=300G

sorted_data=method3.sorted.data.txt.gz
if [ ! -f  $sorted_data ];then
  (
  #for f in $test_ds
  for f in `find . -name "*.gz"`
  do
    id=`echo $f | ruby -ne 'puts $_.match(/(\d+)/)[1]'`
    gzip -cd $f | grep -v -P "^#" | grep -v INDEL | sed 's/Chr//' | awk -v id=$id '{ if ($6>200) print id" "$1" "$2}'
  done
  ) | sort -S$sort_mem -T/space1/tmp -k2,2 -k3,3n | gzip -c > $sorted_data
fi

gzip -cd $sorted_data | go run ./doit.go 4 > results.method3.txt
