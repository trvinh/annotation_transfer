#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Std;

### check the orthology assignment for KO-refProts if they are supported by Inparanoid orthologs
### 21.02.2017

sub usage {
    my $msg = shift;
    print "example: perl compareInpa.pl -h sacce_unknown.list.NEW.KO -i inparanoid/sacce_strict_inpar/sacce.INPAOUT\n";
    print "-h\tKO summary file\n";
    print "-i\tinparanoid output file\n";
    die $msg."\n";
}

# global variables
our($opt_h,$opt_i);
getopts('h:i:');

# sanity checks;
my $hamfasIN = ($opt_h) ? $opt_h : usage("ERROR: No hamstr ortholog list given\n");
my $inpaIN = ($opt_i) ? $opt_i : usage("ERROR: No inparanoid ortholog list given\n");

### get INPARANOID orthologs
open(INPA,$inpaIN) || die "Cannot open $inpaIN!\n";
my %inpa;
foreach my $line(<INPA>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# sacce:100	aspni:7861	neucr:5140	schpo:433
	$inpa{$tmp[0]} = $line;
}
close (INPA);

### parse HamFAS orthologs and compare with INPA orthologs
open(HAM,$hamfasIN) || die "Cannot open $hamfasIN!\n";
open(OUT,">$hamfasIN.noINPA");
open(WITH,">$hamfasIN.withINPA");

my $totalOrtho = 0;
my $noInpaOrtho = 0;
my $totalGroup = 0;	# total groups that have hamstr orthos
my $diffGroup = 0;	# numbers of groups that contains noInpaOrthos
my $groupWoKO = 0;	# numbers of groups have no KO after removing all noInpaOrthos

foreach my $line(<HAM>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# sacce:1227	K13700;0.999838591;nemve:13902	K01259;0.999761533;mge:MG_020	K08726;0.996185722;neucr:8608

	my $seed = shift(@tmp);

	my $flag = 0;	# to mark a group if it has noInpaOrthos
	my @newLine;

	if(scalar(@tmp) > 0){
		# check if this seed protein has any inparanoid ortholog
		if($inpa{$seed}){
			my $inpaOrtho = $inpa{$seed};
			foreach my $prot(@tmp){
				if($prot =~ /;/){
					my @protTMP = split(/;/,$prot);	# K13700;0.999838591;nemve:13902

					# check if this hamfas ortholog NOT supported by inparanoid
					unless($inpaOrtho =~ /$protTMP[2]\t*/){
		#				print "(1) NO INPA FOR $prot with $seed!\n";<>;
						push(@newLine,$prot);
						$noInpaOrtho ++;
						$flag = 1;
					}
					$totalOrtho ++;
				}
			}
		} else {
	#		print "(2) NO INPA FOR @tmp\n";<>;
			my $tmp = join("\t",@tmp);
			push(@newLine,$tmp);

			$noInpaOrtho += scalar(@tmp);
			$totalOrtho += scalar(@tmp);
			$flag = 1;
		}

		$totalGroup ++;
		# count group contains noInpaOrthos
		if($flag == 1){
			$diffGroup ++;
		}

		# count group contains ALL noInpaOrthos
		if(scalar(@newLine) == scalar(@tmp)){
			$groupWoKO ++;
		}
	}

	### print outputs
	my $newLine = join("\t",@newLine);
	$newLine .= "\t";
	# new ortholog groups after removing noInpaOrthos
#	if(scalar(@newLine) < scalar(@tmp)){
		print WITH $seed;
		foreach(@tmp){
			unless($newLine =~ /$_\t/){
				print WITH "\t$_";
			}
		}
		print WITH "\n";
#	}

	# print noInpaOrthos
	if(scalar(@newLine) > 0){
		$newLine =~ s/\t$//;
		print OUT $seed,"\t",$newLine,"\n";
	}

	
}
close (HAM);

print "FINISHED! Output: $hamfasIN.noINPA AND $hamfasIN.withINPA\n";
print "$noInpaOrtho / $totalOrtho (",$noInpaOrtho/$totalOrtho,") orthologs are not supported by Inparanoid\n";
print "$diffGroup / $totalGroup (",$diffGroup/$totalGroup,") are effected\n";
print "$groupWoKO / $totalGroup have no annotation anymore\n";

close (OUT);
close (WITH);

exit;
