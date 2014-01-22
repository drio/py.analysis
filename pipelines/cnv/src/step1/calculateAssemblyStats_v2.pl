#!/usr/bin/env perl
use strict;
use warnings;


# Output files
my $outfile = "stats.txt";
my $logfile = "stats.log";


# Define counters
my %gNtCount = (A => 0, C => 0, G => 0, T => 0, N => 0, l_case => 0);
my $gBp = 0;
my $gBpOther = 0;
my $line = 0;
my $chr;

# Read fasta files
while (<STDIN>) {
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
            print STDERR "$chr at line $line contains a base that is not A, G, T, or N (base = $base)\n";
            next;
        }
        $gBp++;
        $gNtCount{$base}++;
    }
}
print STDOUT "$chr\t$gBp\t".$gNtCount{"A"}."\t".$gNtCount{"C"}."\t".$gNtCount{"G"}."\t".$gNtCount{"T"}."\t".$gNtCount{"N"}."\t".$gNtCount{"l_case"}."\t$gBpOther\n";

