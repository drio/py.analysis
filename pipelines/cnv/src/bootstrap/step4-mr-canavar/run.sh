#!/bin/bash
#
source ../common.sh

sam_dir=`readlink -e ../step3`

if [ ".$1" == "." ]
then
  $step4_sh $fasta_masked $sam_dir | $csv2submit -
else
  echo "test.."
fi
