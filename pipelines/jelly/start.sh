#!/bin/bash

genome=`pwd | ruby -F/ -ane 'puts $F[-2]'`
round=`pwd | ruby -F/ -ane 'puts $F[-1]'`
seed=$genome.$round

cmd="Jelly.py setup ./Protocol.xml 2> logs/setup.log"

echo '$cmd' | submit -s $seed.first | bash > first.jid.txt
cat first.jid.txt | submit -s $seed.mapping -f - "Jelly.py mapping ./Protocol.xml 2> logs/mapping.log" | bash
