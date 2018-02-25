#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
replace old NCBI IDs by new ones for traceability matrix
13.11.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl modify_ncbiID_trace.pl -i sacce_unk.trace\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");

my $newID = "/Users/trvinh/work/OLD/koAnnotation/toolValidation/sacce_rbh_checkco/refSpec.list";
my $oldID = "/Users/trvinh/work/OLD/koAnnotation/toolValidation/sacce_rbh_checkco/refSpec.list.OLD";

open(NEW,$newID);
my %new;
foreach my $line(<NEW>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# Ashbya gossypii (Eremothecium gossypii)	ncbi33169	ago
	$new{$tmp[2]} = $tmp[1];
}
close (NEW);

open(OLD,$oldID);
my %old;
foreach my $line(<OLD>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# aae	ncbi224324
	$old{$tmp[1]} = $tmp[0];
}
close (OLD);

open(IN,$inFile) || die "Cannot open $inFile\n";
my @in = <IN>;
close (IN);
my $in = join("",@in);

my $new = $in;
foreach my $key (keys %old){
	if($new =~ /$key[\t\n]/){
		print "oldID = $key - newID = $new{$old{$key}}\n";
		my $lastChar = substr($&,length($&)-1,1);
		$new =~ s/$&/$new{$old{$key}}$lastChar/;
		unless($new{$old{$key}}){
			print "NO NEW ID for $old{$key}\n";
		}
		unless($old{$key}){
			print "NO OLD ABBR for $key\n";
		}
	} else {
		print "NO $key\n";<>;
	}
}

open(OUT,">$inFile.modID");
print OUT $new;
close (OUT);
exit;
