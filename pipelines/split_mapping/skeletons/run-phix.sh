#!/bin/bash

set -e

[ $(which sapi.py) ] || (echo "sapi.py not found"; exit 1)

sapi.py check_env

curl -L http://cl.ly/1f0q3X2L373U/phix.tar.bz2 | tar -jx

end="-x"
bam="`pwd`/phix/phix.bam"
fa="`pwd`/phix/phix.fa"
#url="-u http://localhost:5000/sapi/api/v1.0/samples"
url=""

sapi.py $url -i FOO -b $bam fastqc $end
sapi.py $url -i FOO -b $bam init -n 40000 $end
sapi.py $url -i FOO -b $bam  -n 40000 splits $end
sapi.py $url -i FOO -b $bam -f $fa sais $end
sapi.py $url -i FOO -b $bam  -f $fa sampe $end
sapi.py $url -i FOO -b $bam -f $fa merge $end
sapi.py $url -i FOO dups $end
sapi.py $url -i FOO stats $end
