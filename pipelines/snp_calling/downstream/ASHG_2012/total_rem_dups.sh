#!/bin/bash
set -e
#set -x

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")
for fc in "raw" "NON_SYNONYMOUS_CODING" "START_GAINED" "START_LOST" "STOP_GAINED" "STOP_LOST" "UTR_5_PRIME"
do
  echo -ne "$fc: "
  cat *.$fc.vcf | awk '{print $1" "$2}' | sort -T/tmp -S8G -k1,1 -k2,2n | uniq | wc -l
done
IFS=$SAVEIFS
