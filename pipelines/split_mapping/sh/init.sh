#!/bin/bash

mkdir -p init
cd init

BAM=$1
NRS=$2

total=`samtools view $BAM |wc -l`
cof=`echo "$total/$NRS" |bc`
rem=`echo "$total%$NRS" |bc`
nsplits=$cof
if [ $rem -ne 0 ];then
    nsplits=$[$nsplits+1]
fi
echo $cof > cof.txt
echo $rem > rem.txt
echo $nsplits > nsplits.txt
