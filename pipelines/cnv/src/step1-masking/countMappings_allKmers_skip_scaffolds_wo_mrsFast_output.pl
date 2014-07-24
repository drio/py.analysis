#!/usr/bin/env perl
# vim: set ts=2 noet sw=2:
use strict;
use warnings;


my $al_file   = $ARGV[0];
my $kmer_file = $ARGV[1];
my $chr       = $ARGV[2];

# Define counters
my %kmersSeqs;
my %kmersCounts;
my $kmers = 0;
my $pair = 0;
my $kmerId;
my $kmerSeq;
my $kmerLen;

my $brokenLines = 0;
my $brokenLinesDiffer = 0;

# Hash the kmer file [id] -> seq; [id] -> 0
open FILE, $kmer_file or die "Could not open $kmer_file: $!";
print STDERR "\n";
print STDERR "# Defining hash with kmers for $chr\n";
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
    print STDERR "***** the number of kmers read is not equal to the number of elements in the hash! *****\n";
}
close FILE;

# Count the number of hits per each kmer
# Check that the alignment matches matches the chrm we are working on
print STDERR "# Counting the number of times each K-mers from $chr is mapped\n";
open FILE, $al_file or die "Could not open $al_file: $!";
while (<FILE>) {
		my $_chr = substr $_, -1;
    if ($_chr ne "\n") {
        print STDERR "skipping broken line: $_\n";
				$brokenLines++;
				next;
		}
    chomp();

    # Expected line format:
    # 11:134658435-134658471  0       1       923942  255     36M     *       0       0       GGTGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAA    *       NM:i:2  MD:Z:10G8A16
		# Only increment the counter (kmer hits) if have a valid alignment record
    my @row = split("\t", $_);
		my $len = @row;
		if ($len == 3) {
			my @subrow = split(":", $row[0]);
			if ($chr eq $subrow[0]) {
				$kmersCounts{$row[0]}++;
			} else {
				$brokenLines++;
				print STDERR "***** chromosomes differ! $chr != $subrow[0] *****\n";
			}
		} else {
			$brokenLines++;
		}
}
print STDERR "# of brokenLines: $brokenLines\n";

# Dump the number of hits per kmer (number of places where the kmer maps to)
foreach (keys %kmersCounts) {
    next if ($kmersCounts{$_} <= 1);
    my @row = split("[:-]", $_);
    print STDOUT "$row[0]\t$row[1]\t$row[2]\t" . $kmersSeqs{$_} . "\t" . $kmersCounts{$_} . "\n";
}
close FILE;
