#!/bin/bash
#

fasta=$1
fasta_anno_masked=$2

[ ".$fasta" == "." ] && echo "Need fasta file." && exit 1
[ ".$fasta_anno_masked" == "." ] && echo "Need name of output fasta file." && exit 1

if [ ! -f $fasta_anno_masked ]
then
    mkdir -p beds
    #curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | gzip -cd - | cut -f 6-8 > beds/rMask.bed
    #curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz | gzip -cd - | cut -f 2-4 > beds/trf.bed
    #curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz | gzip -cd - | awk '{print $2"\t"$3"\t"$4}' > beds/gaps.bed
    cat beds/*.bed | sort -T/space1/tmp -S2g -k1,1 -k2,2n > all.bed
    bedtools maskfasta -fi $fasta -bed all.bed -fo $fasta_anno_masked
    mrsfast --index $fasta_anno_masked
    rm -f all.bed
fi
