#!/bin/bash
#
[ ! -f "./config.sh" ] && echo "I can't find ./config.sh" && exit 1
src_dir=$(dirname $(readlink -f $0))
source ./config.sh
[ ! -f $fasta ] && echo "I can't find $fasta" && exit 1
[ ! -f $fasta_anno_masked ] && echo "I can't find $fasta" && exit 1


pipe="bash" # "" to not submit or "bash" to submit to cluster
th_over="20"


chrm_info="chrm_info.bed"
fasta_anno_masked="fasta_anno_masked.fa"
fa_without_over="anno_kmer_masked.fa"

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
    -e "select chrom, size from $ucsc_genome.chromInfo" | \
    grep chr | awk '{print $1"\t"$2"\t"$2}' | grep -v size > $chrm_info

[ ! -f $chrm_info ] && echo "I can't find $chrm_info" && exit 1
rm -rf deps; mkdir deps


# You may want to run this manually
###################################
cmd="$src_dir/mask_rtf_simple_gaps.sh $fasta $fasta_anno_masked"
submit -e -s "gen_masked_ref" $cmd | $pipe >> deps/1.txt


while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
		cmd="$src_dir/makeIntervalsBED_v2.pl $chrm $len 36 5"
		cat deps/1.txt | submit -f - -e -s "locations.$chrm" $cmd | $pipe >> deps/2.txt
done < $chrm_info


while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    len=$(echo -e $line | cut -f2 -d' ')
    bed="${chrm}_k36_step5_intervals.bed"
    out="${chrm}.reads.fa"
    cmd="fastaFromBed -fi $fasta_anno_masked -bed $bed -fo $out"
    cat deps/2.txt | submit -e -f - -s "extract.$chrm" "$cmd" | $pipe >> deps/3.txt
done < $chrm_info


while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    reads="${chrm}.reads.fa"
    out="${chrm}.sam"
		# FIXME
		# Hack to avoid failing dependency.
		# There is one read file with all N's and mrsfast exit != 0
    cmd="mrsfast --search $fasta_anno_masked --seq $reads -o $out -e 2; exit 0"
    cat deps/3.txt | submit -e -f - -m 16G -s "map_kmers.$chrm" "$cmd" | $pipe >> deps/4.txt
done < $chrm_info


while read line;do
    chrm=$(echo -e $line | cut -f1 -d' ')
    kmer_file="${chrm}.reads.fa"
    alignments="${chrm}.sam"
		# FIXME
		# Hack to avoid failing dependency.
		# There is one read file with all N's and mrsfast exits with value different to 0
    cmd="$src_dir/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl \
			$alignments $kmer_file $chrm > counts.$chrm; exit 0"
    cat deps/4.txt | submit -e -m 16G -f - -s "count.$chrm" "$cmd" | $pipe >> deps/5.txt
done < $chrm_info


cmd="$src_dir/process_counts.sh $th_over"
cat deps/5.txt | submit -e -f - -s "working_on_counts" $cmd | $pipe >> deps/6.txt


cmd="maskFastaFromBed -fi $fasta_anno_masked -bed all.$th_over.counts -fo $fa_without_over"
cat deps/6.txt | submit -e -f - -s "mask_final" "$cmd" | $pipe
