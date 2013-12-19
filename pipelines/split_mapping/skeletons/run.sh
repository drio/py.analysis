#!/bin/bash

[ $(which sapi.py) ] || (echo "sapi.py not found"; exit 1)
[ ! -f config.sh ] && (echo "No config.sh found." ; exit 1);

source config.sh

#clean
#sapi.py -i $id -b $bam init -n $nrps -x
end="-s csv"
(
#sapi.py -i $id -b $bam fastqc $end
[ ! -d "./splits" ] && sapi.py -p /space1/tmp -i ${id} -b $bam  -n $nrps -p $tmp -u $url splits $end
[ ! -d "./sais" ] && sapi.py -p /space1/tmp -i ${id} -t$t -b $bam -f $fa sais -p $tmp -u  $url  $end
[ ! -d "./sampe" ] && sapi.py -p /space1/tmp -i ${id} -b $bam  -f $fa sampe -p $tmp -u $url  $end
[ ! -d "./merge" ] && sapi.py -p /space1/tmp -i ${id} merge -f $fa  -b $bam  -p $tmp -u $url $end
[ ! -d "./dups" ] && sapi.py -p /space1/tmp -i ${id} dups -p $tmp -u $url $end
[ ! -d "./stats" ] && sapi.py -p /space1/tmp -i ${id} stats -p $tmp -u $url $end
) | csv2cluster.py -
