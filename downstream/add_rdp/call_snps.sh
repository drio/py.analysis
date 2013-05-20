#!/bin/bash
#
ref="/stornext/snfs6/rogers/drio_scratch/baboon.diversity/papio_papio_analysis/attempt2/ref/baboon.fa"
bams="/stornext/snfs6/rogers/drio_scratch/baboon.diversity/papio_papio_analysis/attempt2/bams"
#bams="`pwd`/test_bams"
/stornext/snfs3/rogers/drio_scratch/snp-calling-pipe/snpcalling-pipe.sh $bams wgs $ref ""
