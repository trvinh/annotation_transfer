#!/usr/bin/env perl

use WWW::Mechanize;
use strict;
use Cwd;
use Getopt::Std;

### RUN fetch KO information from list of genome IDs
### 02.09.2016

sub usage {
    my $msg = shift;
    print "example: perl run_fetchKO.pl -i 4koAnnotation.list\n";
    print "-i\tinput genome IDs list\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

my $input = ($opt_i) ? $opt_i : usage("ERROR: No list of genome IDs given!\n");

open(IN,$input) || die "Cannot open $input!\n";
my @in = <IN>;
close (IN);

#my @nodeList = ("leonard","leslie","raj","beverly","will","alex","hawking","amy","penny","stuart","bernadette");
#my $c = 1;
foreach my $line(@in){
	chomp($line);
	my @tmp = split(/\t/,$line);

	my $cmd = "perl fetchKO_fromKEGG_new.pl -i ".$tmp[1]." -o ".$tmp[0];
#	system($cmd);

#	my $i = int($c/2);
	my $cluster_run = "echo \"$cmd\" | qsub -V -S /bin/bash -cwd -j y -N $tmp[0] -r n -q workstation.q\@leonard.izn-ffm.intern";
	system($cluster_run);
#	$c ++;
}

exit;
