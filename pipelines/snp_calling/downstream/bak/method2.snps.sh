#!/bin/bash


#Chr1    2749    .       T       C       222     .       DP=31;VDB=0.0399;AF1=1;AC1=2;DP4=0,0,18,11;MQ=59;FQ=-114        GT:PL:GQ        1/1:255,87,0:99
#Chr1    2764    .       G       A       222     .       DP=28;VDB=0.0398;AF1=1;AC1=2;DP4=0,0,17,8;MQ=59;FQ=-102 GT:PL:GQ        1/1:255,75,0:99
#Chr1    2933    .       C       A       222     .       DP=24;VDB=0.0345;AF1=1;AC1=2;DP4=0,0,11,10;MQ=60;FQ=-90 GT:PL:GQ        1/1:255,63,0:99
GS="2900000000"
#  echo -ne "$i: "; ns=`zcat $i | grep -v "#"| grep -v INDEL | wc -l`; echo "scale=10; ($ns/$GS)*1000"| bc

sort_mem=400G
find . -name "*.gz" | xargs -i zcat {} | grep -v "#" | grep -v INDEL | awk '{print $1" "$2}' | sed 's/Chr//' | sort -S$sort_mem -T/space1/tmp -k1,1 -k2,2n |\
uniq -c > pos_n_samples.method2.txt
