#!/bin/perl
use strict;
use warnings;


# Input file
my $infile = $ARGV[0];
my $assembly = $ARGV[1];
my $version = $ARGV[2];
my $chrId = $ARGV[3];
open FILE, $infile or die "Could not open $infile: $!";


# Output files
my $outfile = "/scratch/primate/Assemblies/rhesusMac2/CNV/$version/AssemblyStats/tmp/tmp.$chrId.assemblyStats.txt";
my $logfile = "/scratch/primate/Assemblies/rhesusMac2/CNV/$version/AssemblyStats/assemblyStats.log";
open OUTFILE, ">>$outfile" or die;
open LOGFILE, ">>$logfile" or die;
#print OUTFILE "chromosome\ttotal\tA\tC\tG\tT\tN\tother\n";
#print "chromosome\ttotal\tA\tC\tG\tT\tN\tother\n";

# Directory
#my @FILES = glob("/scratch/primate/Assemblies/CanFam2/chromFaMasked_kmerMasked/*.fa.kmer.masked");
#foreach (@FILES) {
	
	# Input file
	#my $infile = $_;


	# Define counters
	my %gNtCount = (A => 0, C => 0, G => 0, T => 0, N => 0, l_case => 0);
	my $gBp = 0;
	my $gBpOther = 0;
	my $line = 0;
	my $chr;

	# Read fasta files
	while (<FILE>) {
		$line++;
		# Exclude mitocondrial
		#next if ($_ =~ /^>chrM/);
		# Get chromosome name	
		if ($_ =~ /^>(.+)/) {
			$chr = $1;
			print "Reading chromosome $chr\n";
			next;
		}
		# Add to counters
		for (my $i = 0; $i < (length($_) - 1); $i++) {
			my $base = substr($_, $i, 1);
			# Count lower-case bases
			if ($base =~ /[acgtn]/) {
				$gNtCount{"l_case"}++;
				$gBp++;
				next;
			}
			# Count bases		
			unless ($base =~ /[ACGTN]/) {
				$gBpOther++;
				print LOGFILE "$chr at line $line contains a base that is not A, G, T, or N (base = $base)\n";
				next;
			}
			$gBp++;
			$gNtCount{$base}++;				
		}
	}
	print OUTFILE "$chr\t$gBp\t".$gNtCount{"A"}."\t".$gNtCount{"C"}."\t".$gNtCount{"G"}."\t".$gNtCount{"T"}."\t".$gNtCount{"N"}."\t".$gNtCount{"l_case"}."\t$gBpOther\n";
	#print "$chr\t$gBp\t".$gNtCount{"A"}."\t".$gNtCount{"C"}."\t".$gNtCount{"G"}."\t".$gNtCount{"T"}."\t".$gNtCount{"N"}."\t$gBpOther\n";
#}


close OUTFILE;
close LOGFILE;



