#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my %opt=();

my ($program, $pversion, $pdescription, $pusage);

$program = $0;
$program =~ s/^.*\///;

$pversion='1';
$pdescription = "$program (ver:$pversion)  ";
$pusage="$0 ";

if (! defined $ARGV[0]){
  print "
    DESCRIPTION $pdescription 
    USAGE $pusage 
      -f1 input file 1 with coordinates chr st e
      -h1 there is a header in -f1
	  optional, by default it is not
      -f2 optional, input file 2 with coordinates chr st e
      -h2 there is a header in -f2
	  optional, by default it is not
      -d  optional, by default 0, maximum distance to merge blocks if they are separated by a length smaller than d
      -outF output file \n";
  exit;
}

if (! &GetOptions(\%opt, "f1:s","f2:s","d:i","h1","h2", "outF:s")){
  die "  Command line unsucessfully parsed!\n";
}


#CHECKING THE COMMAND and INPUTS
$opt{'f1'} ||= die "  Insert a file with the first coordinates to merge, with the command -f1\n";
#$opt{'f2'} ||= die "  Insert a file with the second coordinates to merge, with the command -f2\n";
$opt{'d'} ||=0;
$opt{'outF'} ||= die "  Insert the output file, with the command -outF\n";


my %cds;
open (IN, $opt{'f1'}) or die "could not open $opt{'f1'}\n";
if (defined $opt{'h1'}){
	<IN>;
}
while (my $line = <IN>){
	chomp($line);

	next if ($line eq "");
	my @cols= split (/\t/,$line);
	if (!defined $cds{$cols[0]}){
		$cds{$cols[0]}=[];
	}
	if ($cols[2]<$cols[1]){
		my $aux = $cols[1];
		$cols[1] = $cols[2];
		$cols[2] = $aux;
	}
	push @{$cds{$cols[0]}}, [$cols[1],$cols[2]];
}
close IN;

if (defined $opt{'f2'}){
	if (defined $opt{'h2'}){
		<IN>;
	}
	open (IN, $opt{'f2'}) or die "could not open $opt{'f2'}\n";
	while (my $line = <IN>){
		chomp($line);

		next if ($line eq "");
		my @cols= split (/\t/,$line);
		if (!defined $cds{$cols[0]}){
			$cds{$cols[0]}=[];
		}
		if ($cols[2]<$cols[1]){
			my $aux = $cols[1];
			$cols[1] = $cols[2];
			$cols[2] = $aux;
		}
		push @{$cds{$cols[0]}}, [$cols[1],$cols[2]];
	}
	close IN;
}

open (OUT, ">".$opt{'outF'}) or die "could not open $opt{'outF'}\n";
foreach my $q (sort(keys %cds)){
	@{$cds{$q}} = sort{ $a->[0] <=> $b->[0]} @{$cds{$q}};
                                                
	my $numInterv = scalar(@{$cds{$q}});
	foreach my $i (0..($numInterv-2)){
		if ($cds{$q}->[$i+1]->[0] < 1+$cds{$q}->[$i]->[1]+$opt{'d'}){
			$cds{$q}->[$i+1]->[0] = $cds{$q}->[$i]->[0];
			if($cds{$q}->[$i+1]->[1] < $cds{$q}->[$i]->[1]){
				$cds{$q}->[$i+1]->[1] = $cds{$q}->[$i]->[1];
			}                       
			$cds{$q}->[$i] = undef;
  		}                       
	}
	foreach my $elem (@{$cds{$q}}){
		next unless defined $elem;

		print OUT "$q\t$elem->[0]\t$elem->[1]\n";
	}

}
close OUT;


