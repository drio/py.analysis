#!/bin/perl
use strict;
use warnings;


#print "***** Check assembly paths! *****\n";
#exit;


# Get parammeters from script arguments
my $cName = $ARGV[0];
my $cLen = $ARGV[1];
my $assembly = $ARGV[2];



# Size of kmers and step
my $k = 36;
my $step = 5;


print "# defining $k -mers(step = $step) for chromosome $cName (length = $cLen)\n";
my $intervals = "/scratch/primate/Assemblies/rhesusMac2/CNV/intervals/" . $cName . "_k" . $k . "_step" . $step . "_intervals.bed";
open INTERVALS, ">$intervals" or die;
my $b = 1;
my $e = $b + $k - 1;
while ($e <= $cLen) {
	print INTERVALS "$cName\t" . ($b-1) . "\t$e\n";
	$b = $b + $step;
	$e = $b + $k - 1;
}	
close INTERVALS;	
