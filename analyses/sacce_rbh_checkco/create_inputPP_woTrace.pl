#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
create input file for phyloprofile
13.11.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl create_inputPP.pl -i seqID_hamfasOnly.list -f sacce_anno.list.NEW.KO -r refSpec.list\n";
		print "-i list of sequence IDs\n";
		print "-f file contains KO annotations\n";
		print "-r list of abbr names and IDs of reference species\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f,$opt_r);
getopts('i:f:r:');

my $listIn = ($opt_i) ? $opt_i : usage("ERROR: No input list given\n");
my $koIn = ($opt_f) ? $opt_f : usage("ERROR: No input KO annotation given\n");
my $refIn = ($opt_r) ? $opt_r : usage("ERROR: No reference species file given\n");

### get species info
open(REF,$refIn) || die "Cannot open $refIn\n";
my %refSpec;	# $refSpec{ncbiID} = abbrName;
foreach my $line(<REF>){
	chomp($line);
	unless($line =~ /^fullName/){
		my @tmp = split(/\t/,$line);	# Aspergillus nidulans	ncbi162425	aspni
		$refSpec{$tmp[1]} = $tmp[2];
	}
}
close (REF);

### read KO annotation file
open(KO,$koIn) || die "Cannot open $koIn!\n";
my %ko;
foreach my $line(<KO>){
	chomp($line);
	my @tmp = split(/\t/,$line);
	$ko{$tmp[0]} = $line;
}
close (KO);

### MAIN
open(IN,$listIn) || die "Cannot open $listIn!\n";
open(OUT,">$listIn.phyloprofile");
print OUT "geneID\tncbiID\torthoID\tFAS\n";
# open(MISS,">$listIn.missingTrace");

foreach my $line(<IN>){
  chomp($line);
	my $geneID = $line;

  if($ko{$line}){
    my @info = split(/\t/,$ko{$line});
		my @allSpec = keys %refSpec;
    for(my $i=1; $i<scalar(@allSpec); $i++){
			# print $firstLine[$i]," - ",$refSpec{$firstLine[$i]}," - ";
			my $ncbiID = $allSpec[$i];
			my $orthoID = "NA";
      my $fas = "NA";
			if($ko{$geneID} =~ /\K\d{5};[01]\.\d+;$refSpec{$allSpec[$i]}\:(.)+/){	# K01655;0.99893297801;ago:AGOS_ADR107W
				my $hit = $&;
				my @hitTMP = split(/;/,$hit);
				$orthoID = $hitTMP[@hitTMP-1];
				$fas = $hitTMP[1];
			}
			# print "$orthoID#$fas\n";
			unless($orthoID eq "NA"){
				print OUT "$geneID\t$ncbiID\t$orthoID\t$fas\n";
			}
			# <>;
    }
  } else {
    print "$line missing\n";
		# print OUT "$geneID\t$ncbiID\t$orthoID\t$fas\t$traceInfo[$i]\n";
  }
}
close (IN);
close (OUT);

exit;
