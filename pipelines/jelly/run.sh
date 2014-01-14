#!/bin/bash
set -e
mkdir -p logs
for i in setup      mapping   support extraction assembly  output
do
  echo "Jelly.py $i ./Protocol.xml 2> logs/$i.log"
done

