#!/usr/bin/env perl
use strict;
use warnings;

#  run "$SRC_DIR/makeIntervalsBED_v2.pl $chrm $len $kmer_size $step"

# Get parammeters from script arguments
my $cName = $ARGV[0];
my $cLen = $ARGV[1];
my $k = $ARGV[2];
my $step = $ARGV[3];

print "# defining $k -mers(step = $step) for chromosome $cName (length = $cLen)\n";
my $intervals = "./" . $cName . "_k" . $k . "_step" . $step . "_intervals.bed";
open INTERVALS, ">$intervals" or die;
my $b = 1;
my $e = $b + $k - 1;
while ($e <= $cLen) {
	print INTERVALS "$cName\t" . ($b-1) . "\t$e\n";
	$b = $b + $step;
	$e = $b + $k - 1;
}
close INTERVALS;
