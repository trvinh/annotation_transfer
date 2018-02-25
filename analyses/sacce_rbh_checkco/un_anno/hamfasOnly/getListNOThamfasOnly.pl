#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
get list of proteins that are annotated not only by hamfas
08.11.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl analyzeLengthPfamOrig.pl -i sacce_unknown.list.NEW.KO -h sacce_unknown.list.NEW.KO_hamfasOnly\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_h);
getopts('i:h:');

my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $hamfasOnlyIn = ($opt_h) ? $opt_h : usage("ERROR: No hamfas_only file given\n");

open(HAM,$hamfasOnlyIn) || die "Cannot open $hamfasOnlyIn!\n";
my @ham = <HAM>;
close (HAM);
my $ham = join("",@ham);

open(IN,$inFile) || die "Cannot open $inFile!\n";
open(OUT,">$inFile.others");

foreach my $line(<IN>){
	chomp($line);
	if($line =~ /\t/){
		my @tmp = split(/\t/,$line);
		unless($ham =~ /$tmp[0]\t/){
			print OUT $line,"\n";
		}
	}
}

close (IN);
close (OUT);

exit;
