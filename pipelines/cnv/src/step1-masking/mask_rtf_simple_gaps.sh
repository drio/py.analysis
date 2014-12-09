#!/bin/bash
#

fasta=$1
fasta_anno_masked=$2
ucsc_code=$3
src_dir=$(dirname $(readlink -f $0))
mrsfast="$src_dir/../../bin/mrsfast"

[ ".$fasta" == "." ] && echo "Need fasta file." && exit 1
[ ".$fasta_anno_masked" == "." ] && echo "Need name of output fasta file." && exit 1
[ ".$ucsc_code" == "." ] && echo "Need ucsc code." && exit 1

if [ ! -f $fasta_anno_masked ]
then
    mkdir -p beds
    [ ! -f beds/rMask.bed ] && curl http://hgdownload.soe.ucsc.edu/goldenPath/${ucsc_code}/database/rmsk.txt.gz | gzip -cd - | cut -f 6-8 > beds/rMask.bed
    [ ! -f beds/trf.bed ] && curl http://hgdownload.soe.ucsc.edu/goldenPath/${ucsc_code}/database/simpleRepeat.txt.gz | gzip -cd - | cut -f 2-4 > beds/trf.bed
    [ ! -f beds/gaps.bed ] && curl http://hgdownload.soe.ucsc.edu/goldenPath/${ucsc_code}/database/gap.txt.gz | gzip -cd - | awk '{print $1"\t"$2"\t"$3}' > beds/gaps.bed
    [ ! -f all.bed ] && cat beds/*.bed | sort -T/space1/tmp -S2g -k1,1 -k2,2n > all.bed
    [ ! -f $fasta_anno_masked ] && bedtools maskfasta -fi $fasta -bed all.bed -fo $fasta_anno_masked
    $mrsfast --index $fasta_anno_masked
fi
