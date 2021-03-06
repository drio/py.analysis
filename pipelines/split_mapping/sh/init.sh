#!/bin/bash


mkdir -p init
cd init

BAM=$1
NRS=$2
URL=$3
ID=$4

[ ! -f ${BAM}.bai ] && [ ! -w `dirname ${BAM}` ] && echo "Bam dir not writable. I can't create index." && exit 1
[ ! -f ${BAM}.bai ] && samtools index $BAM

date > time.txt
total=`samtools idxstats $BAM | awk '{s+=$3+$4} END {print s}'`
cof=`echo "$total/$NRS" |bc`
rem=`echo "$total%$NRS" |bc`
nsplits=$cof
if [ $rem -ne 0 ];then
    nsplits=$[$nsplits+1]
fi
echo $cof > cof.txt
echo $rem > rem.txt
echo $nsplits > nsplits.txt

signal.py $URL $ID init

date >> time.txt
