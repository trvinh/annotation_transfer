#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
get list of length, number of pfam domains and origin for
list of annotated proteins (e.g. sacce_unknown.list.NEW.KO)
08.11.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl analyzeLengthPfamOrig.pl -i sacce_unknown.list.NEW.KO_hamfasOnly -a FACT_annotation/koAnnotation/sacce/pfam.xml\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_a);
getopts('i:a:');

my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $pfamIN = ($opt_a) ? $opt_a : usage("ERROR: No PFAM annotation file given\n");

### get age of reference species where KO annotation come from
my %age;
$age{'sacce'}=1;
$age{'ago'}=1;
$age{'aspni'}=1;
$age{'canal'}=1;
$age{'neucr'}=1;
$age{'schpo'}=1;

$age{'rno'}=2;
$age{'hsa'}=2;
$age{'mmu'}=2;

$age{'cel'}=3;
$age{'dme'}=3;
$age{'dre'}=3;
$age{'nemve'}=3;
$age{'monbr'}=3;
$age{'cho'}=3;
$age{'pfa'}=3;
$age{'trybr'}=3;
$age{'ath'}=3;
$age{'ehi'}=3;

$age{'ape'}=4;
$age{'mja'}=4;

$age{'aae'}=5;
$age{'bsu'}=5;
$age{'eco'}=5;
$age{'hpy'}=5;
$age{'lla'}=5;
$age{'mge'}=5;
$age{'mtu'}=5;
$age{'nme'}=5;
$age{'syn'}=5;

my %type = (
	"1" => "fungi",
	"2" => "mammals",
	"3" => "eukaryotes",
	"4" => "archaea",
	"5" => "bacteria",
	"6" => "others"
);

### PFAM annotation file
open(ANNO,$pfamIN) || die "Cannot open $pfamIN!\n";
my @anno = <ANNO>;
close (ANNO);
my $anno = join("",@anno);
my @allAnno = split(/protein>/,$anno);

my %len;
my %domain;
foreach my $annoProt (@allAnno){
	my @tmp = split(/\n/,$annoProt);
	my $id = "";
	for(my $i=0; $i < scalar(@tmp); $i++){
		if($tmp[$i] =~ /id=(.)+/){
			my $hit = $&;
			my @hitTMP = split(/\s/,$hit);
			$id = $hitTMP[0];	$id =~ s/id=//; $id =~ s/\"//g;
			my $len = $hitTMP[1];	$len =~ s/length=//; $len =~ s/\"//g; $len =~ s/>//;
			# print $id,"\t",$len;<>;
			$len{$id} = $len;
			last;
		}
	}
	my @c = $annoProt =~ /feature type/g;
	my $count = @c;
	# print $count;<>;
	$domain{$id} = $count;
}

### get info for input protein list
open(IN,$inFile) || die "Cannot open $inFile!\n";
open(OUT,">$inFile.LengthDomain");
print OUT "seqID\tlength\tdomains\torigin\n";

foreach my $line(<IN>){
	chomp($line);

	# get seqID
	my @tmp = split(/\t/,$line);		# sacce:3479		K12388;0.96312050102;rno:83576
	my $seqID = shift(@tmp);
	print OUT $seqID,"\t";

	# print seq length
	if($len{$seqID}){
		print OUT $len{$seqID},"\t";
	} else {
		print OUT "noLen\t";
	}
	# print number of Pfam domains
	print OUT $domain{$seqID},"\t";

	# get and print origin of orthologs
	my $minAge = 6;
	foreach my $item(@tmp){
		if($item =~ /K/){
			if($seqID eq "sacce:3471"){
				print $item,"\t";<>;
			}
			# print $item,"\n";
			my @itemTMP = split(/;/,$item);
			my @refSpec = split(/:/,$itemTMP[@itemTMP-1]);
			# print "AGE:$age{$refSpec[0]}\n";<>;
			if($minAge > $age{$refSpec[0]}){
				$minAge = $age{$refSpec[0]};
				if($seqID eq "sacce:3471"){
					print "$minAge - $refSpec[0]\n";
				}
			}
		}
	}
	print OUT $type{$minAge},"\n";
}
close (IN);
close (OUT);
print "FINISHED! Check output $inFile.LengthDomain\n";
exit;
