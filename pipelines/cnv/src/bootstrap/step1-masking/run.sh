#!/bin/bash
#

source ../common.sh

$step1_sh $fasta | $csv2submit -
