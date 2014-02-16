#!/usr/bin/env perl
use strict;
use warnings;


my $kmer_file = $ARGV[0];
my $chr       = $ARGV[1];

# Define counters
my %kmersSeqs;
my %kmersCounts;
my $kmers = 0;
my $pair = 0;
my $kmerId;
my $kmerSeq;
my $kmerLen;

# Open file
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

# Open file
print STDERR "# Counting the number of times each K-mers from $chr is mapped\n";
while (<STDIN>) {
    chomp();
    my @row = split("\t", $_);
    $kmersCounts{$row[0]}++;
    my @subrow = split(":", $row[0]);
    unless ($chr eq $subrow[0]) {
        print STDERR "***** chromosomes differ! $chr != $subrow[0] *****\n";
    }
}

# Print to output file
foreach (keys %kmersCounts) {
    next if ($kmersCounts{$_} <= 1);
    my @row = split("[:-]", $_);
    print STDOUT "$row[0]\t$row[1]\t$row[2]\t" . $kmersSeqs{$_} . "\t" . $kmersCounts{$_} . "\n";
}
