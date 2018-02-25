#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
get list of proteins from KO file
05.12.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl getLines.pl -i ids.list -f sacce_unknown.list.KO\n";
    print "-i\tinput ID list\n";
    print "-f\toriginal file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f);
getopts('i:f:');

my $inList = ($opt_i) ? $opt_i : usage("ERROR: No input ID list given\n");
my $inFile = ($opt_f) ? $opt_f : usage("ERROR: No original file given\n");

### read and parse original file
open(ID,"$inList") || die "Cannot open $inList!\n";
my @ids = <ID>;
close (ID);
my $ids = join("",@ids);
$ids .= "\n";

### search those IDs in original file
open(OUT,">$inFile.filtered");
open(IN,$inFile) || die "Cannot open $inFile!\n";
foreach my $line(<IN>){
	chomp($line);
	my @tmp = split(/\t/,$line);
	if($ids =~ /$tmp[0]\n/){
		print OUT $line,"\n";
	}
}
close (IN);
close (OUT);

print "FINISHED! Check output $inFile.filtered\n";
exit;
