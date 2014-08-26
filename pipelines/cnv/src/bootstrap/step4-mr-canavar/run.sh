#!/bin/bash
#
source ../common.sh

sam_dir=`readlink -e ../step3-mapping`

if [ ! -f $fasta_masked_pad ];then
    echo "Cannot find padding masked ref"
    exit
fi

if [ ".$1" == "." ]
then
  $step4_sh $fasta_masked_pad $sam_dir | $csv2submit -
else
  echo "test.."
fi
