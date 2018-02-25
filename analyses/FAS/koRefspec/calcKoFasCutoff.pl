#!/usr/bin/env perl
use warnings;
use strict;
use Cwd;
use Cwd 'abs_path';
use Getopt::Std;
use POSIX;
use List::Util qw(sum);

=desc
calculate KO threshold
02.09.2016
=cut

sub usage {
    my $msg = shift;
    print "example: perl cutoff.pl -i koRefspec.list.FAS\n";
    print "-i\tCalculated FAS score list of all KO refspec proteins\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');
my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");

open(IN,$inFile) || die "Cannot open $inFile!\n";
my @in = <IN>;
close (IN);

my %koFas;
print "Parsing input file...";
foreach my $line(@in){
	chomp($line);
#	print $line;<>;
	my @tmp = split(/\t/,$line);	# K00001	ath:AT5G24760	ath:AT1G64710	0.99831437044
	unless($koFas{$tmp[0]}){
		$koFas{$tmp[0]} = $tmp[@tmp-1];
	} else {
		$koFas{$tmp[0]} .= ";".$tmp[@tmp-1];
	}
}
print "done!\n";

### calculate median scores
open(MEAN,">koThreshold_mean.list");
open(MEDIAN,">koThreshold_median.list");

print "KO\tmean\tmedian\n";
my $c = 0;
foreach my $ko(sort keys %koFas){
	my @allFas = split(/;/,$koFas{$ko});
	my $mean = mean(@allFas);
	my $median = median(@allFas);
#	print $ko,"\t",$mean,"\t",$median,"\n";<>;
	print MEAN $ko,"\t",$mean,"\n";
	print MEDIAN $ko,"\t",$median,"\n";

	$c ++;
	print $c," / ",scalar(keys %koFas),"\n";
}

print "FINISHED!! Check outputs in\nkoThreshold_mean.list and koThreshold_median.list\n";
exit;

sub median{
	return sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

sub mean {
	return sum(@_)/@_;
}
