#!/bin/bash
#
# vim: set ts=2 noet ft=sh sw=2:
# CNV pipeline: First step
# Masking over-represented parts of the genome
#
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BN=`basename $0`

for c in index intervals extract map count mask hist help;do
	if [ ".$RE_COMMANDS" == "." ];then
		RE_COMMANDS="^$c\$"
	else
		RE_COMMANDS="$RE_COMMANDS|^$c\$"
	fi
done

error() {
	local msg=$1
	echo $1 >&2
	exit 1
}

usage() {
	[ ".$1" != "." ] && echo "ERROR: $1"
	cat <<EOF
Usage: $BN <command> [options]

Valid commands (use $BN help <command> for command details):

	index    : index reference genome.
	intervals: generate genomic intervals for the kmers.
	extract  : extract kmers from reference.
	map      : map kmers to genome.
	count    : count number of hits per kmer.
	mask     : mask genome using over represented kmers coordinates.
	hist     : make histogram of counts.
	help     : get command help
EOF
	exit
}

parse_args() {
	[ $# -lt 1 ] && usage "Invalid number of arguments."
	command=$1
	shift

  [[ $command =~ $RE_COMMANDS ]] || usage "Invalid command."
	[[ $command == "help" ]] && help_for=$1

	kmer_size=36
	step=5
	debug=0

	while getopts a:b:c:l:f:k:s:i:r:m:t:d flag; do
		[ "$OPTARG" == "-" ] && OPTARG="/dev/stdin"
		case $flag in
			'a') alignments=$OPTARG;;
			'b') int_bed=$OPTARG;;
			'c') chrm=$OPTARG;;
			'l') len=$OPTARG;;
			'f') fasta_genome=$OPTARG;;
			'k') kmer_size=$OPTARG;;
			's') step=$OPTARG;;
			'i') genome_info=$OPTARG;;
			'r') reads=$OPTARG;;
			'm') kmer_file=$OPTARG;;
			't') counts=$OPTARG;;
			'd') debug=1;;
			?) help;;
		esac
	done
}

check_env() {
	for b in R grep cut awk mrsfast maskFastaFromBed fastaFromBed
	do
		`which $b &>/dev/null` || error "$b not in path."
	done
}

check_var() {
	[[ $# -ne 2 || ".$1" == "." ]] && usage $2
}

check_file() {
	[[ ! -f $1 ]] && usage "Cannot find file: $1"
}

help() {
	if [[ $help_for =~ $RE_COMMANDS ]];then
		help_$help_for
	else
		usage "That command does not exist."
		exit 0
	fi
}

run() {
	[ $debug == 1 ] && echo $1 || $1
}

help_defaults() {
	cat <<EOF

kmer_size: Size of the kmers [$kmer_size]
step_size: Skip S bp for continous kmer [$step]
EOF
}

help_genome_info() {
	cat <<EOF

Format of genome info file (bed):
chrm end end
chr1 45600 45600
...

You may find this useful:
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from rheMac2.chromInfo" | \
grep chr | awk '{print \$1"\t"\$2"\t"\$2}' | grep -v size > genome_info.bed
EOF
}

help_index() {
	echo "Usage: $BN index -f <fasta.fa>"
}

index() {
	check_var $fasta_genome "Need fasta genome."
	run "mrsfast --index $fasta_genome --ws 12 ${fasta_genome}.index"
}

help_intervals() {
	cat <<EOF
Usage: $BN intervals -c <chrm> -l <length> [-k <kmer_size> -s <step_size>]
EOF
	help_defaults
}

intervals() {
	check_var $chrm "Need chrm name"
	check_var $len "Need lenght of the chrm"
  run "$SRC_DIR/makeIntervalsBED_v2.pl $chrm $len $kmer_size $step"
}

help_extract() {
	echo "Usage: $BN extract -f <genome_fasta_file> -b <intervals>"
}

extract() {
	check_var $fasta_genome "Need fasta genome"
	check_var $int_bed "Need kmer intervals"
	run "fastaFromBed -fi $fasta_genome -bed $int_bed -fo stdout"
}

help_map() {
	echo "Usage: $BN map -f <genome_fasta_file> -r <reads_file>"
}

map() {
	check_var $fasta_genome "Need fasta genome"
	check_var $reads "Need input reads"
	run "mrsfast --search $fasta_genome --seq $reads -o /dev/stdout -e 2"
}

help_count() {
	echo "Usage: $BN counts -m <kmer_file> -a <alignments> -c <chrm>"
}

count() {
	check_var $kmer_file "Need kmer file."
	check_var $alignments "Need alignments."
	check_var $chrm "Need chrm."
	run "${SRC_DIR}/countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl $alignments $kmer_file $chrm"
}

help_mask() {
	echo "Usage: $BN mask -f <fasta_genome> -t <over_counts_bed>"
}

mask() {
	check_var $fasta_genome "Need fasta genome file."
	check_var $counts "Need counts file."
	run "maskFastaFromBed -fi $fasta_genome -bed $counts -fo /dev/stdout"
}

help_hist() {
	echo "Usage: $BN hist -t <over_counts_bed>"
}

hist() {
	check_var $counts "Need counts file."
	[ $counts == "/dev/stdin" ] && counts=stdin
	debug=1 # FIXME: this command doesn't execute properly when used with run()
	        # That's why I force debugging
	run "xvfb-run R CMD BATCH --no-restore --no-save '--args $counts' ${SRC_DIR}/makeHistogramCounts.R"
}

###################
# Main
###################
check_env
parse_args $*
$command
