#!/bin/bash
#
source ../common.sh

input_dir=`readlink -e ../step2-read-prep`

if [ ".$1" == "." ]
then
  $step3_sh $input_dir $fasta_masked | $csv2submit -
else
  echo "test.."
  exit 0
fi
