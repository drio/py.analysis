#!/bin/perl
use strict;
use warnings;


# Get parammeters from script arguments
my $cName = $ARGV[0];
# $ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from rheMac2.chromInfo"
my $cLen = $ARGV[1];

# Size of kmers and step
my $k = 36;
my $step = 5;

print STDERR "# defining $k -mers(step = $step) for chromosome $cName (length = $cLen)\n";
my $intervals = $cName . "_k" . $k . "_step" . $step . "_intervals.bed";
my $b = 1;
my $e = $b + $k - 1;
while ($e <= $cLen) {
	print STDOUT "$cName\t" . ($b-1) . "\t$e\n";
	$b = $b + $step;
	$e = $b + $k - 1;
}
