#!/bin/bash

[ $(which sapi.py) ] || (echo "sapi.py not found"; exit 1)

curl -L http://cl.ly/1f0q3X2L373U/phix.tar.bz2 | tar -jx

end="-x"
bam="`pwd`/phix/phix.bam"
fa="`pwd`/phix/phix.fa"

sapi.py -i FOO -b $bam fastqc $end
sapi.py -i FOO -b $bam init -n 40000 $end
sapi.py -i FOO -b $bam  -n 40000 splits $end
sapi.py -i FOO -b $bam -f $fa sais $end
sapi.py -i FOO -b $bam  -f $fa sampe $end
sapi.py -i FOO -b $bam -f $fa merge $end
sapi.py -i FOO dups $end
sapi.py -i FOO stats $end
