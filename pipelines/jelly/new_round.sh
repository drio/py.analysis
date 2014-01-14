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
mkdir -p $dir/reference
mkdir -p $dir/logs

cp -r ../boilerplate/* ./$dir/
rm -f $dir/new_round.sh
ln -s `pwd`/$prev/jelly.out.fasta ./$dir/reference/input.fasta
rm -f $dir/reference/ml.fasta
rm -rf ./$dir/reads
ln -s `pwd`/$prev/reads ./$dir/

cat ./$dir/Protocol.xml | \
  sed 's:__REF__:$ref:g' | \
  sed 's:__OUT__:$base:g' | \
  sed 's:__BASE__:$base:g' > _tmp
mv _tmp ./$dir/Protocol.xml

EOF
