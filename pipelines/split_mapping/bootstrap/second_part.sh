#!/bin/sh

source config.sh

#sapi.py -i $id -b $bam init -n $nrps -x
end="-s csv"
(
#sapi.py -i $id -b $bam fastqc $end
[ ! -d "./splits" ] && sapi.py -i ${id}_${seed} -b $bam  -n $nrps splits $end
[ ! -d "./sais" ] && sapi.py -i ${id}_${seed} -t$t -b $bam -f $fa sais $end
[ ! -d "./sampe" ] && sapi.py -i ${id}_${seed} -b $bam  -f $fa sampe $end
[ ! -d "./merge" ] && sapi.py -i ${id}_${seed} merge -f $fa  -b $bam $end
[ ! -d "./dups" ] && sapi.py -i ${id}_${seed} dups $end
[ ! -d "./stats" ] && sapi.py -i ${id}_${seed} stats $end
) | csv2cluster.py -
