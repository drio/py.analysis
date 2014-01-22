#!/bin/perl
use strict;
use warnings;


my $assembly = "rheMac2";


# Input files
my $kmersFile = "/scratch/primate/Assemblies/rhesusMac2/CNV/kmerCountsWithMrsFast_more20placements.txt";
my @mrsFastOutFiles = glob "/scratch/primate/Assemblies/rhesusMac2/CNV/mrsFast/*_k36_step5_intervals.map";
my @fastaFiles = glob "/scratch/primate/Assemblies/rhesusMac2/CNV/chromFaMasked/chr*";


# Output files
my $bedFile = "/scratch/primate/Assemblies/rhesusMac2/CNV/intervalsToBeMasked.bed";
open BED, ">$bedFile" or die;


# Define a hash with over-represented kmers
print "\n# Reading K-mers placed >20 times from $kmersFile\n";
open KMERS, $kmersFile or die "Could not open $kmersFile: $!";
my %kmers;
#my $i = 0;
while (<KMERS>) {
	#$i++;
	#if ($i % 1e5 == 0) {
	#	print "$i\n";
	#}
	my @row = split("\t", $_);
	my $kmerId = "$row[0]:$row[1]-$row[2]";
	$kmers{$kmerId} = $row[3];
}


# Get placements for the kmers
print "\n# Get placements for K-mers in:\n";
foreach (@mrsFastOutFiles) {
	# Skip *.fai files
	next if ($_ =~ m/\.fai$/);
	print "\t$_\n";
	#next;
	open MRSFASTOUT, $_ or die "Could not open $_: $!";
	while (<MRSFASTOUT>) {
		my @row = split("\t", $_);
		# Select only K-mers that defined previously as over-placed
		if (exists $kmers{$row[0]}) {
			# Check that sequences match
			# Consider sequences in complementary/reverse orientation
			my $seq;
			if ($row[1] == 0) {
				$seq = $row[9];
			}			
			elsif ($row[1] == 16) {
				$seq = reverse $row[9];
				$seq =~ tr/ACGT/TGCA/;
			}
			unless ($kmers{$row[0]} eq $seq) {
				print "***** different sequences! *****\n";
			}
			# Print to BED
			# Note that BED format uses 0-start format
			print BED "$row[2]\t" . ($row[3]-1) . "\t" . ($row[3]-1+36) . "\n";
			# Add to the hash the information on whether the K-mer is placed
		}
	}
}

close BED;
