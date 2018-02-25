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
    print "example: perl create_inputPP.pl -i seqID_hamfasOnly.list -f sacce_anno.list.NEW.KO -t sacce_anno_hamfasOnly.trace -r refSpec.list\n";
		print "-i list of sequence IDs\n";
		print "-f file contains KO annotations\n";
		print "-t computed traceability matrix\n";
		print "-r list of abbr names and IDs of reference species\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f,$opt_t,$opt_r);
getopts('i:f:t:r:');

my $listIn = ($opt_i) ? $opt_i : usage("ERROR: No input list given\n");
my $traceIn = ($opt_t) ? $opt_t : usage("ERROR: No input trace matrix given\n");
my $koIn = ($opt_f) ? $opt_f : usage("ERROR: No input KO annotation given\n");
my $refIn = ($opt_r) ? $opt_r : usage("ERROR: No reference species file given\n");

### get species info
open(REF,$refIn) || die "Cannot open $refIn\n";
my %refSpec;	# $refSpec{ncbiID} = abbrName;
foreach my $line(<REF>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# Aspergillus nidulans	ncbi162425	aspni
	$refSpec{$tmp[1]} = $tmp[2];
}
close (REF);

### get trace info
open(TRACE,$traceIn) || die "Cannot open $traceIn!\n";
my @trace = <TRACE>;
close (TRACE);
my $firstLine = shift(@trace);
chomp($firstLine);
my @firstLine = split(/\t/,$firstLine);
my %trace;
foreach my $line(@trace){
  chomp($line);
  my @tmp = split(/\t/,$line);
  $trace{$tmp[0]} = $line;
}

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
print OUT "geneID\tncbiID\torthoID\tFAS\ttraceability\n";
open(MISS,">$listIn.missingTrace");

foreach my $line(<IN>){
  chomp($line);
  if($trace{$line}){
    my @traceInfo = split(/\t/,$trace{$line});
    my $geneID = $traceInfo[0];
    for(my $i=1; $i<scalar(@traceInfo); $i++){
			# print $firstLine[$i]," - ",$refSpec{$firstLine[$i]}," - ";
			my $ncbiID = $firstLine[$i];
			my $orthoID = "NA";
      my $fas = "NA";
			if($ko{$geneID} =~ /\K\d{5};[01]\.\d+;$refSpec{$firstLine[$i]}\:(.)+/){	# K01655;0.99893297801;ago:AGOS_ADR107W
				my $hit = $&;
				my @hitTMP = split(/;/,$hit);
				$orthoID = $hitTMP[@hitTMP-1];
				$fas = $hitTMP[1];
			}
			# print "$orthoID#$fas\n";
			print OUT "$geneID\t$ncbiID\t$orthoID\t$fas\t$traceInfo[$i]\n";
			# <>;
    }
  } else {
    print MISS "$line\n";
  }
}
close (IN);
close (OUT);
close (MISS);

exit;
