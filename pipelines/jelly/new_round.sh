#!/bin/bash

curr=`ls -d round* | tail -1 | sed 's/round//g'`
if [ $curr == "." ];then
  echo "Not in a jelly round dir?"
  exit 1
fi
prev=round$curr

dir=round$[$curr+1]
base=`pwd`/$dir
ref=`pwd`/$dir/reference/input.fasta

cat << EOF
mkdir -p $dir/logs

rm -rf ./$dir/reads
ln -s `pwd`/$prev/reads ./$dir/

mkdir -p $dir/reference
ln -s `pwd`/$prev/jelly.out.fasta ./$dir/reference/input.fasta

cp -r ../boilerplate/stats ./$dir/

cat `pwd`/$prev/Protocol.xml | sed 's:$prev:$dir:g' > ./$dir/Protocol.xml
EOF
