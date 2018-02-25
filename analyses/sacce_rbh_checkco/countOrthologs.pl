#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
count number of orthologs
07.12.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl countOrthologs.pl -i sacce_unk.list.NEW.KO\n";
    print "-i\tinput KO annotated file\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input annotated file given\n");

open(OUT,">$inFile.orthoCount");
open(IN,$inFile) || die "Cannot open $inFile!\n";
foreach my $line(<IN>){
	if($line =~ /K\d{5}/){
		chomp($line);
		my @tmp = split(/\t/,$line);
		print OUT $tmp[0],"\t",scalar(@tmp)-1,"\n";
	}
}
close (IN);
close (OUT);
print "FINISHED! Check output $inFile.orthoCount\n";
exit;
