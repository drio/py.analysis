#!/bin/bash

set -e
for i in setup      mapping   support extraction assembly  output
do
  if [ ! -f "logs/$i.log" ];then
    ./run.sh | grep $i
    exit
  fi
done

echo "It seems you are done?"

