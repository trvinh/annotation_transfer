#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;
use POSIX;
use Scalar::Util qw(looks_like_number);

=desc
compare annotated KOs with original annotation
23.08.2016
=cut

sub usage {
    my $msg = shift;
    print "example: perl compareKoalaHamfas.pl -i hamfas.KO -a /home/vinh/Desktop/data/project/KEGG_annotation/koMap\n";
    print "-i: input KO annotation file\n";
    print "-a: original KO annotation folder\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_a);
getopts('i:a:');

# sanity checks;
my $fileIN = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $koPath = ($opt_a) ? $opt_a : usage("ERROR: No KO annotation folder given\n");

###### MAIN ######
open(IN,$fileIN) || die "Cannot open $fileIN!!\n";
my @input = <IN>;
close (IN);

my %koOrig;	# original KO
my $parsedKO = "";	# list of parsed KO species

my @totalRelevant;
my @totalRetrieved;
my @relevantRetrieved;
my @irrelevantRetrieved;
my @missing;

open(IRRE,">$fileIN.irrelevant");
open(RELFAS,">$fileIN.relevantFAS");
open(IRREFAS,">$fileIN.irrelevantFAS");
my $relfas_tag = 0;
my $irrelfas_tag = 0;

foreach my $line(@input){
	chomp($line);

	# get species name/ID
	my $flag = 0;	# $flag = 1: input file is NOT hamfas result
	my $specID = "";
	if($line =~ /\|/){
		$flag = 1;		# ANNO_10001|sacce|sacce:2188	K00888
		my @lineTMP = split(/\t/,$line);
		my @specTMP = split(/\|/,$lineTMP[0]);
		$specID = $specTMP[1];
	} else {
		my @lineTMP = split(/\:/,$line);		# sacce:2522      K03245;0.99719658654;dme:Dmel_CG1213    K08139;0.99299695971;ago:AGOS_AFL205C	
		$specID = $lineTMP[0];
	}

	# open original annotation file if not yet parsed
	unless($parsedKO =~ /$specID/){
		if(-e "$koPath/$specID.ko"){
			open(KO,"$koPath/$specID.ko");
			my @koOrig = <KO>;
			close (KO);
			foreach my $koline(@koOrig){
				chomp($koline);
				my @tmp = split(/\t/,$koline);	# ape:APE_0166.1  K00798  2.5.1.17
				$koOrig{$tmp[0]} = $tmp[1];
			}
			$parsedKO .= ";".$specID;
		} else {
			print "NO original annotation for $specID\n";<>;
		}
	}

	### COMPARE RESULT WITH ORIGINAL ANNOTATION
	my $origKO = "NONE";

	# if there is any KO predicted
	if($line =~ /\tK\d{5}/){
		push(@totalRetrieved,$line);

		my @lineTMP = split(/\t/,$line);		# sacce:2522      K03245;0.99719658654;dme:Dmel_CG1213    K08139;0.99299695971;ago:AGOS_AFL205C
		# get protein ID
		my $protID = shift(@lineTMP);
		if($flag == 1){		# ANNO_10001|sacce|sacce:2188	K00888
			my @protTMP = split(/\|/,$protID);
			$protID = pop(@protTMP);
		}
		# if this protein is already annotated 
		if($koOrig{$protID}){
			push(@totalRelevant,$protID);

			# check if this is a relevant retrieved (predicted KO = original KO)
			$origKO = $koOrig{$protID};
			my $result = $protID."\t".$origKO."\t";
			$result .= join("\t",@lineTMP);
			if($line =~ /$koOrig{$protID}/){
				push(@relevantRetrieved,$result);
				# get FAS
				if($line =~ /$koOrig{$protID}\;(.)+?\;/){
					my $fasTMP = $&;
					my @fasTMP = split(/\;/,$fasTMP);
					print RELFAS $fasTMP[1],"\n";
					$relfas_tag = 1;
				}
			} else {
				print IRRE "DIFFERENT for $protID: old $koOrig{$protID} - new $line\n";
				push(@irrelevantRetrieved,$result);
				# get FAS
				if($lineTMP[0] =~ /\;/){
					my @fasTMP = split(/\;/,$lineTMP[0]);
					print IRREFAS $fasTMP[1],"\n";
					$irrelfas_tag = 1;
				}
			}
		} else {
			print IRRE "NEW KO for $protID: $line\n";
			my $result = $protID."\t".$origKO."\t";
			$result .= join("\t",@lineTMP);
			push(@irrelevantRetrieved,$result);
			# get FAS
			if($lineTMP[0] =~ /\;/){
				my @fasTMP = split(/\;/,$lineTMP[0]);
				print IRREFAS $fasTMP[1],"\n";
				$irrelfas_tag = 1;
			}
		}
	} else {
#		print $line;<>;
		my $protID = $line;
		$protID =~ s/\t*$//;
		if($flag == 1){			
			my @protTMP = split(/\|/,$line);
			$protID = pop(@protTMP);
		}
		if($koOrig{$protID}){
			my $result = $protID."\t".$koOrig{$protID};
			push(@totalRelevant,$protID);
			push(@missing,$result);
		}
	}
}
close (IRRE);

my $recall = scalar(@relevantRetrieved)/scalar(@totalRelevant);
my $precision = scalar(@relevantRetrieved)/scalar(@totalRetrieved);
my $fscore = 2*$recall*$precision/($recall+$precision);
print "\ntotalRelevant: ",scalar(@totalRelevant),"\n";
print "total retrieved: ",scalar(@totalRetrieved),"\n";
print "relevant retrieved: ",scalar(@relevantRetrieved),"\n";
print "irrelevant retrieved: ",scalar(@irrelevantRetrieved),"\n";
print "missing: ",scalar(@missing),"\n\n";
print "recall: $recall\nprecison: $precision\nF-score: $fscore\n\n";

if($irrelfas_tag == 0){system("rm $fileIN.irrelevantFAS");}
if($relfas_tag == 0){system("rm $fileIN.relevantFAS");}

exit;
