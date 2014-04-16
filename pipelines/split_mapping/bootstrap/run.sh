#!/bin/sh

(submit -s _one 'source config.sh; sapi.py -i $id -b $bam -n $nrps init -x' | bash) | \
        submit -s _two -f - "./second_part.sh | bash" | bash
