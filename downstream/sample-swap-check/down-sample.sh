#!/bin/bash
#
# vim: set ts=2 sw=2 noet paste:
#
# shuffle a fastq stream (stdin)
#
n_of_records=10000
[ ".$1" != "." ] && n_of_records=$1


# paste f1.fastq f2.fastq
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | \
	shuf  | \
	head -$n_of_records | \
	sed 's/\t\t/\n/g'
# awk '{print $1 > "file1.fastq"; print $2 > "file2.fatsq"}'
