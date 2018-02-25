#!/usr/bin/env perl
use strict;
use warnings;

### create all pairwise proteins for koRefspec list

open(IN,"koRefspec.list");
my @in = <IN>;
close (IN);

open(OUT,">koRefspec.list.pairs");
open(SINGLE,">koRefspec.list.single");
foreach my $line(@in){
	chomp($line);
	if(length $line > 2){
		my @lineTMP = split(/\t/,$line);	# @lineTMP = koID prot1 prot2 prot3
		my $koID = shift(@lineTMP);

		if(scalar(@lineTMP) < 2){
			print SINGLE "$line\n";
		} else {
			foreach my $prot1(@lineTMP){
				foreach my $prot2(@lineTMP){
					unless ($prot1 eq $prot2){
						print OUT $koID,"\t",$prot2,"\t",$prot1,"\n";
					}
				}
			}
		}
	}
}

close (OUT);
print "Finished!! Output files are koRefspec.list.pairs & koRefspec.list.single\n";

exit;
