#!/bin/perl
use strict;
use warnings;


my $assembly = "rheMac2";


# Define paths
my @kmersFiles = glob "/scratch/primate/Assemblies/rhesusMac2/CNV/intervalsFasta/chr*_k36_step5_intervals.fa.masked";
my @mappingsFiles = glob "/scratch/primate/Assemblies/rhesusMac2/CNV/mrsFast/chr*_k36_step5_intervals.map";

# Output file
my $out = "/scratch/primate/Assemblies/rhesusMac2/CNV/kmerCountsWithMrsFast.txt";
open OUT, ">$out" or die;


# Define those chromosomes for which mrsFAST output mappings
my %mapped;
foreach (@mappingsFiles) 
{
	my $chr = $_;
	$chr =~ s/_k36_step5_intervals.map//;
	$chr =~ s/\/scratch\/primate\/Assemblies\/rhesusMac2\/CNV\/mrsFast\///;

	unless (exists $mapped{$chr}) 
	{
		$mapped{$chr} = 1;
	}
}
	
# Set hashes
foreach (@kmersFiles) 
{
	# Get chromosome name
	my $chr = $_;
	$chr =~ s/_k36_step5_intervals.fa.masked//;
	$chr =~ s/\/scratch\/primate\/Assemblies\/rhesusMac2\/CNV\/intervalsFasta\///;

	unless (exists $mapped{$chr}) 
	{
		next;
	}

	# Define counters	
	my %kmersSeqs;
	my %kmersCounts;
	my $kmers = 0;
	my $pair = 0;	
	my $kmerId;
	my $kmerSeq;
	my $kmerLen;
	# Open file
	open FILE, $_ or die "Could not open $_: $!";
	print "\n";
	print "# Defining hash with kmers for $chr\n";
	while (<FILE>) {
		# Get kmer identifier
		if ($_ =~ /^>(.*)/) {$kmerId = $1; $pair++;}
		# Get kmer sequence
		else {$kmerSeq = $_; chomp($kmerSeq); $kmerLen = length $kmerSeq; $pair++;}
		# Add to hashes
		if ($pair == 2) {
			$pair = 0;
			# Skip if all kmer are N's
			next if ($kmerSeq =~ /N{$kmerLen}/);
			$kmersSeqs{$kmerId} = $kmerSeq;
			$kmersCounts{$kmerId} = 0;
			$kmers++;				
		}
	}
	# Check that the number of kmers read = number of elements in the hash
	unless ($kmers == scalar keys %kmersSeqs) 
	{
		print "***** the number of kmers read is not equal to the number of elements in the hash! *****\n";
	}
	# Count K-mers
	my $mapFile = "/scratch/primate/Assemblies/rhesusMac2/CNV/mrsFast/" . $chr . "_k36_step5_intervals.map";
	# Open file
	open FILE, $mapFile or die "Could not open $mapFile: $!";
	print "# Counting the number of times each K-mers from $chr is mapped\n";
	while (<FILE>) {
		my @row = split("\t", $_);
		$kmersCounts{$row[0]}++;
		my @subrow = split(":", $row[0]);
		unless ($chr eq $subrow[0]) {
			print "***** chromosomes differ! *****\n";
		}
	}
	# Print to output file
	foreach (keys %kmersCounts) {
		next if ($kmersCounts{$_} <= 1);
		my @row = split("[:-]", $_);
		print OUT "$row[0]\t$row[1]\t$row[2]\t" . $kmersSeqs{$_} . "\t" . $kmersCounts{$_} . "\n";
	}
}

close OUT;
