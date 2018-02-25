#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;
use Cwd 'abs_path';

=desc
create ref ko list
K0001	prot1	prot2	prot3
K0008	prot4	prot5
...
v1.0 (08.08.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl createKoRefList.pl -i refSpec.list -d KEGG_annotation/ko\n";
    print "-i List of ref species\n";
    print "-d Folder contains KOs of those refSpec\n";
    die $msg."\n";
}

### check input
our($opt_i,$opt_d);
getopts('i:d:');
my $listIN = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $koFol = ($opt_d) ? $opt_d : usage("ERROR: No input folder given\n");

### open input file
open(IN,$listIN) || die "Cannot open $listIN\n";
my @in = <IN>;
close (IN);

my $koPath = abs_path($koFol);

### MAIN
my %koList;	# $koList{$koID} = prot1  prot2   prot3...

my $c = 1;
foreach my $spec(@in){
	chomp($spec);
	if(length $spec > 2){
		unless(-e "$koPath/$spec.ko"){
			print "No KO file for $spec in $koPath!\n";<>;
		} else {
			open(KO,"$koPath/$spec.ko");
			my @koFile = <KO>;
			close (KO);

			### read this KO file
			foreach my $line(@koFile){
				chomp($line);
				if(length $line > 2){
					my @lineTMP = split(/\t/,$line);	# dre:678543      K16330  3.2.-.-;2.7.1.83
	
					### save each KO identifier together with its proteins into %koList
					unless($koList{$lineTMP[1]}){
						$koList{$lineTMP[1]} = $lineTMP[0];
					} else {
						$koList{$lineTMP[1]} .= "\t".$lineTMP[0];
					}
				}
			}

		}
		print $c," - ",$spec,"\n"; $c++;
	}
}

### print output
open(OUT,">koRefspec.list");
foreach (sort keys %koList){
	print OUT $_,"\t",$koList{$_},"\n";
}
close (OUT);

print "Finished!!! Output file is koRefspec.list\n";
exit;
