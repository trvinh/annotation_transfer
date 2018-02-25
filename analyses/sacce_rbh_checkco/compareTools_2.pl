#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
analyse the results of unknown protein annotation
1) calculate the mean FAS score of new annotated KO
2) compare between HamFAS, blastKOALA and KAAS
01.09.2016
print annotated KO (if available)
=cut

sub usage {
    my $msg = shift;
    print "example: perl compareTool.pl -h hamfas.KO -k kaas.KO -b blastKOALA.KO -a koMap\n";
    print "-h\tKO annotation obtained from HamFAS\n";
    print "-k\tKO annotation obtained from KAAS\n";
    print "-b\tKO annotation obtained from blastKOALA\n";
    print "-a\tFolder contains original KO lists\n";
    die $msg."\n";
}

# global variables
our($opt_h,$opt_k,$opt_b,$opt_a);
getopts('h:k:b:a:');

my $hamIN = ($opt_h) ? $opt_h : usage("ERROR: No HamFAS annotation given\n");
my $kaasIN = ($opt_k) ? $opt_k : usage("ERROR: No KAAS annotation given\n");
my $koalaIN = ($opt_b) ? $opt_b : usage("ERROR: No blastKOALA annotation given\n");
my $koPath = ($opt_a) ? $opt_a : usage("ERROR: No original KO folder given\n");

### variables
my @hamOnly;
my @kaasOnly;
my @koalaOnly;

my @hamKaas;
my @hamKoala;
my @kaasKoala;

my @all;
my @noAnno;

### original KOs
my %koOrig;	# original KOs
my $parsedKO = "";	# list of parsed KO species

### get kaas and blastkoala annotations
my %kaas = getKOList($kaasIN);
my %koala = getKOList($koalaIN);

### HamFAS 
open(HAM,$hamIN) || die "Cannot open $hamIN!\n";
my @ham = <HAM>;
close (HAM);

open(FASOUT,">$hamIN.fasList");
open(HAMONLY,">$hamIN.fasList_hamfasOnly");

foreach my $hamLine(@ham){
	chomp($hamLine);
	my @tmp = split(/\t/,$hamLine);	# sacce:791       K03068;0.97990263401;ampqu_4652_1:I1GE66        K01586;0.0;aae:aq_1208
	my $protID = shift(@tmp);
	my $hamKO = join("#",@tmp);

	# open original annotation file if not yet parsed
	my @protIDtmp = split(/\:/,$protID);		# sacce:2522      K03245;0.99719658654;dme:Dmel_CG1213    K08139;0.99299695971;ago:AGOS_AFL205C	
	my $specID = $protIDtmp[0];

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

	# get original KO if any
	my $oriKO = "";
	if($koOrig{$protID}){$oriKO = $koOrig{$protID};}

	# compare HamFAS to kaas and blastKoala
	if($hamLine =~ /K\d{5}/){
		# get list of all (best) FAS scores, used for FAS score statistical analysis
		my @fasTMP = split(/;/,$tmp[0]);	### best KO annotation for this protein
		if($fasTMP[1] >= 0){
			print FASOUT $fasTMP[1],"\n";
		}

		# check if this proteins has kaas AND blastkoala annotation
		if($kaas{$protID} and $koala{$protID}){
			my $result = $protID."\t".$oriKO."\t".$kaas{$protID}."\t".$koala{$protID}."\t".$hamKO;

			## compare with kaas and blastkoala
			if ($hamLine =~ /$kaas{$protID}/ or $hamLine =~ /$koala{$protID}/){
				if($tmp[0] =~ /$kaas{$protID}/ or $tmp[0] =~ /$koala{$protID}/){
					$result .= "\t1";
				} else {
					$result .= "\t0.5";
				}
			} else {
				$result .= "\t0";
			}

			## compare with original ko
			if($oriKO eq ""){
				$result .= "\tn";
			} else {
				if ($hamLine =~ /$oriKO/){
					if($tmp[0] =~ /$oriKO/){
						$result .= "\t1";
					} else {
						$result .= "\t0.5";
					}
				} else {
					$result .= "\t0";
				}
			}

			push(@all,$result);
		} else {
			# check if this proteins has kaas (BUT NOT blastkoala) annotation
			unless($koala{$protID}){
				if($kaas{$protID}){
					my $result = $protID."\t".$oriKO."\t".$kaas{$protID}."\t".$hamKO;
					if ($hamLine =~ /$kaas{$protID}/){
						if($tmp[0] =~ /$kaas{$protID}/){
							$result .= "\t1";
						} else {
							$result .= "\t0.5";
						}
					} else {
						$result .= "\t0";
					}

					## compare with original ko
					if($oriKO eq ""){
						$result .= "\tn";
					} else {
						if ($hamLine =~ /$oriKO/){
							if($tmp[0] =~ /$oriKO/){
								$result .= "\t1";
							} else {
								$result .= "\t0.5";
							}
						} else {
							$result .= "\t0";
						}
					}

					push(@hamKaas,$result);
				} else {
					my $result = $protID."\t".$oriKO."\t".$hamKO;

					## compare with original ko
					if($oriKO eq ""){
						$result .= "\tn";
					} else {
						if ($hamKO =~ /$oriKO/){
							$result .= "\t1";
						} else {
							$result .= "\t0";
						}
					}

					push(@hamOnly,$result);
					print HAMONLY $fasTMP[1],"\n";
				}
			} else {
				# check if this proteins has NO kaas BUT blastkoala annotation
				unless($kaas{$protID}){
					my $result = $protID."\t".$oriKO."\t".$koala{$protID}."\t".$hamKO;
					if ($hamLine =~ /$koala{$protID}/){
						if($tmp[0] =~ /$koala{$protID}/){
							$result .= "\t1";
						} else {
							$result .= "\t0.5";
						}
					} else {
						$result .= "\t0";
					}

					## compare with original ko
					if($oriKO eq ""){
						$result .= "\tn";
					} else {
						if ($hamLine =~ /$oriKO/){
							if($tmp[0] =~ /$oriKO/){
								$result .= "\t1";
							} else {
								$result .= "\t0.5";
							}
						} else {
							$result .= "\t0";
						}
					}

					push(@hamKoala,$result);
				}
			}
		}
	} else {
		# check if this proteins has kaas AND blastkoala annotation
		if($kaas{$protID} and $koala{$protID}){
			my $result = $protID."\t".$oriKO."\t".$kaas{$protID}."\t".$koala{$protID};
			if ($kaas{$protID} eq $koala{$protID}){
				$result .= "\t1";
			} else {
				$result .= "\t0";
			}

			## compare with original ko
			if($oriKO eq ""){
				$result .= "\tn";
			} else {
				if ($koala{$protID} eq $oriKO){
					$result .= "\t1";
				} else {
					$result .= "\t0";
				}
			}

			push(@kaasKoala,$result)
		} else {
			if($kaas{$protID}){
				my $result = $protID."\t".$oriKO."\t".$kaas{$protID};

				## compare with original ko
				if($oriKO eq ""){
					$result .= "\tn";
				} else {
					if ($kaas{$protID} eq $oriKO){
						$result .= "\t1";
					} else {
						$result .= "\t0";
					}
				}

				push(@kaasOnly,$result);
			} elsif ($koala{$protID}){
				my $result = $protID."\t".$oriKO."\t".$koala{$protID};

				## compare with original ko
				if($oriKO eq ""){
					$result .= "\tn";
				} else {
					if ($koala{$protID} eq $oriKO){
						$result .= "\t1";
					} else {
						$result .= "\t0";
					}
				}

				push(@koalaOnly,$result);
			} else {
				my $result = $protID."\t".$oriKO;
				push(@noAnno,$result);
			}
		}
	}
}
close (FASOUT);
close (HAMONLY);

open(OUT,">stat.out") || die "Cannot create stat.out !\n";

print OUT "###############################\nALL\t",scalar(@all),"\n";
print OUT "KAAS & HamFAS\t",scalar(@hamKaas),"\n";
print OUT "KOALA & HamFAS\t",scalar(@hamKoala),"\n";
print OUT "KAAS & KOALA\t",scalar(@kaasKoala),"\n";
print OUT "HamFAS only\t",scalar(@hamOnly),"\n";
print OUT "KAAS only\t",scalar(@kaasOnly),"\n";
print OUT "BlastKOALA only\t",scalar(@koalaOnly),"\n";
print OUT "No annotation\t",scalar(@noAnno),"\n";
print OUT "###############################\n\n";

print OUT "###### ALL\t",scalar(@all),"\n";
my $all = join("\n",@all);
print OUT $all,"\n";

print OUT "###### KAAS & HamFAS\t",scalar(@hamKaas),"\n";
my $hamKaas = join("\n",@hamKaas);
print OUT $hamKaas,"\n";

print OUT "###### KOALA & HamFAS\t",scalar(@hamKoala),"\n";
my $hamKoala = join("\n",@hamKoala);
print OUT $hamKoala,"\n";

print OUT "###### KAAS & KOALA\t",scalar(@kaasKoala),"\n";
my $kaasKoala = join("\n",@kaasKoala);
print OUT $kaasKoala,"\n";

print OUT "###### HamFAS only\t",scalar(@hamOnly),"\n";
my $hamOnly = join("\n",@hamOnly);
print OUT $hamOnly,"\n";

print OUT "###### KAAS only\t",scalar(@kaasOnly),"\n";
my $kaasOnly = join("\n",@kaasOnly);
print OUT $kaasOnly,"\n";

print OUT "###### BlastKOALA only\t",scalar(@koalaOnly),"\n";
my $koalaOnly = join("\n",@koalaOnly);
print OUT $koalaOnly,"\n";

print OUT "###### No annotation\t",scalar(@noAnno),"\n";
my $noAnno = join("\n",@noAnno);
print OUT $noAnno,"\n";


exit;

sub getKOList{
	my $fileIN = $_[0];
	open(IN,$fileIN) || die "Cannot open $fileIN!\n";
	my @file = <IN>;
	close (IN);

	my %out;
	foreach my $line(@file){
		chomp($line);
		if($line =~ /K\d{5}/){
			my @tmp = split(/\t/,$line);	# ANNO_10009|sacce|sacce:2827     K14007

			my $id = $tmp[0];	# ANNO_10009|sacce|sacce:2827
			if($id =~ /\|/){
				my @idTMP = split(/\|/,$id);
				$id = pop(@idTMP);
			}

			$out{$id} = $tmp[1];
		}
	}

	return %out;
}
