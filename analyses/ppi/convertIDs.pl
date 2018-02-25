#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
convert refseq IDs into inhouse IDs
30.11.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl convertIDs.pl\n";
    die $msg."\n";
}

my $koPath = "/home/vinh/Desktop/data/project/KEGG_annotation/ko";
my $refseqPath = "/home/vinh/Desktop/data/project/KEGG_annotation/refseq";
my $uniprotPath = "/home/vinh/Desktop/data/project/KEGG_annotation/uniprot";
my $faPath = "/share/project/vinh/dataset/allSpeciesFasta";

### get list of species which have KO annotation
my $koSpecList = "/share/project/vinh/dataset/4koAnnotation.list";
open(refList, $koSpecList) || die "Cannot open $koSpecList!\n";
my @refList = <refList>;
close (refList);

#my $allmicrosID = "enche_5516_1;encin_5517_1;enccu_2934_1;nosce_4242_1;entbi_4241_1;vitco_50505;annca_339;antlo_2712;edhae_4124_1;vavcu_5255_1;nempa_5256_1";
#my @refList = split(/;/,$allmicrosID);

my $outPath = "/home/vinh/Desktop/data/project/KEGG_annotation/koMap2";
open(ERROR,">$outPath/error.log");
print ERROR "#PROT_ID\tXREF_ID (not found in) FILE\n";

foreach my $specID (@refList){
	chomp($specID);
	if(length $specID > 1){
		print "###### ",$specID,"\n";
		## open fasta file
		open(FASTA,"$faPath/$specID.fa") || die "Cannot open $faPath/$specID.fa!\n";
		my @faFile = <FASTA>;
		close (FASTA);

		## output file
		unless(-e "$outPath/$specID.ko"){
			open(OUT,">$outPath/$specID.ko");

			## read fasta file
			foreach my $line(@faFile){
				if($line =~ /^>/){
					##### get sequence ID
					chomp(my $protID = $line); $protID =~ s/>$specID\://;
	#				print $protID,"\n";

					##### get KO for this sequence
					my $koID = ""; my $ecID = "";

					if(-e "$koPath/$specID.ko"){
						### get KO
						open(KO,"$koPath/$specID.ko") || die "Cannot open $koPath/$specID.ko !!\n";
						my @koList = <KO>;
						close (KO);
						my $koList = "\n"; $koList = join("",@koList); $koList .= "\n";
						my $keggID = "";

						### for nemve_2309, sacce_2336, aspni_2095, pucgr_1921, neucr_1906 and schpo_2340_1 -> need refseq files
						if($specID =~ /nemve/ || $specID =~ /sacce/ || $specID =~ /aspni/ || $specID =~ /pucgr/ || $specID =~ /schpo/ || $specID =~ /neucr/ || $specID =~ /trybr/){
							# open refseq file
							open(XREF,"$refseqPath/$specID.refseq") || die "Cannot find/open $refseqPath/$specID.refseq !!\n";
							my @xref = <XREF>;
							close (XREF);
							my $xref = "\n";
							$xref .= join("",@xref);
							$xref .= "\n";

							# search for line in xref that contains this ortholog
							if($xref =~ /\n(.)*\t$protID\b/){
	#							print $&,"\n";
								# get XREF_ID
								my @hit = split(/\t/,$&);
								my @xrefIDTMP = split(/\:/,$hit[0]);
								chomp(my $xrefID = $xrefIDTMP[1]);
	#							print $xrefID,"\n";

								# KO for schpo_2340_1 and neucr_1906
								if($specID =~ /schpo/ || $specID =~ /neucr/){
									## open uniprot file
									open(UNI,"$uniprotPath/$specID.uniprot") || die "Cannot find/open $uniprotPath/$specID.uniprot !!\n";
									my @uni = <UNI>;
									close (UNI);
									my $uni = "\n"; $uni .= join("",@uni); $uni .= "\n";
									# edit searching pattern for neucr (from Q7SI03-8337 to Q7SI03)
									if($specID =~ /neucr/){
										my @xrefID_tmp = split(/-/,$xrefID);
										$xrefID = shift(@xrefID_tmp);
									}

									if($uni =~ /$xrefID\b(.)*?\n/){
										my @uniHit = split(/\t/,$&);
										# get kegg ID
										my @uniHit_tmp = split(/\s/,$uniHit[4]);
										$keggID = $uniHit_tmp[@uniHit_tmp-1];
										($koID,$ecID) = getKO($keggID,$koList);
									} else {
										print "Cannot find hit in $specID.uniprot for $xrefID !!\n";
										print ERROR "$protID\t$xrefID\t$specID.uniprot\n";
									#	unless($lackUniProt{$specID}){$lackUniProt{$specID} = $xrefID;}
									#	else {$lackUniProt{$specID} .= "\t".$xrefID;}
									}
								}
								# KO for nemve_2309_1
								elsif($specID =~ /nemve/){
									# parse keggID
									my @xref_tmp = split(/@/,$xrefID);
									$keggID = "NEMVE_v1g".$xref_tmp[0];
									# search for KO (and EC)
									($koID,$ecID) = getKO($keggID,$koList);
	#								if(length $koID > 2) {print "here: $keggID - $koID";<>;}
								}
								# KO for sacce_2336, aspni_2095, pucgr_1921
								elsif($specID =~ /sacce/ || $specID =~ /aspni/ || $specID =~ /pucgr/ || $specID =~ /trybr/){
									# search for KO (and EC)
									($koID,$ecID) = getKO($xrefID,$koList);
	#								if(length $koID > 2) {print "here: $koID";<>;}
								}

							} else {
	#							print "NO XREF HIT!!\n";
							}
						}

						# for lacbi_2940_1, ampqu_4652_1 and monbr_1569_1 (new fasta -> new ID), ngr, psoj
						elsif ($specID =~ /lacbi/){
							my @keggID_tmp = split(/_/,$protID);
							shift(@keggID_tmp);
							$keggID = join("_",@keggID_tmp);
							($koID,$ecID) = getKO($keggID,$koList);
						} elsif($specID =~ /ampqu/){
							my @keggID_tmp = split(/_/,$protID);
							shift(@keggID_tmp);
							$keggID = join("",@keggID_tmp);
							$keggID =~ s/LOC//;
							($koID,$ecID) = getKO($keggID,$koList);
						} elsif($specID =~ /monbr/){
							my @keggID_tmp = split(/_/,$protID);
							shift(@keggID_tmp);
							$keggID = join("",@keggID_tmp);
							$keggID = "MONBRDRAFT_".$keggID;
							($koID,$ecID) = getKO($keggID,$koList);
						} elsif($specID =~ /ngr/){
							my @keggID_tmp = split(/_/,$protID);
							shift(@keggID_tmp);
							$keggID = join("",@keggID_tmp);
							$keggID = "NAEGRDRAFT_".$keggID;
							($koID,$ecID) = getKO($keggID,$koList);
						} elsif($specID =~ /psoj/){
							my @keggID_tmp = split(/_/,$protID);
							shift(@keggID_tmp);
							$keggID = join("",@keggID_tmp);
							$keggID = "PHYSODRAFT_".$keggID;
							($koID,$ecID) = getKO($keggID,$koList);
						}

						### for canal & cre
						elsif($specID =~ /canal/ || $specID =~ /cre/){
							# open uniprot file
							open(UNI,"$uniprotPath/$specID.uniprot") || die "Cannot find/open $uniprotPath/$specID.uniprot !!\n";
							my @uni = <UNI>;
							close (UNI);
							my $uni = "\n"; $uni .= join("",@uni); $uni .= "\n";
#	if($protID eq "Q5A5F2"){
							# search for protID
							my %ko; my %ec;
							if($uni =~ /$protID\b(.)*?\n/){
								my @uniHit = split(/\t/,$&);
								# get kegg ID
								my @uniHit_tmp = split(/\s/,$uniHit[4]);
								my $prefix = "Ca";	# used to search for KEGG gene ID in uniprot file
								if($specID =~ /cre/){ $prefix = "CHLREDRAFT";}
#		print $&,"\n";
								foreach(@uniHit_tmp){
									if($_ =~ /$prefix/){
										my $keggIDtmp = $_;
#		print "HERE: $keggIDtmp\n";
										if($keggIDtmp =~ /\//){
											my @keggIDtmp = split(/\//,$keggIDtmp);
											my @pre = split(/\./,$keggIDtmp);

											foreach my $keggIDtmp2(@keggIDtmp){
												unless($keggIDtmp2 =~ /$pre[0]\b/){
													$keggIDtmp2 = $pre[0].".".$keggIDtmp2;
												}

												unless($keggID =~ /$keggIDtmp2\;/){$keggID .= $keggIDtmp2.";";}
												my ($koID_tmp,$ecID_tmp) = getKO($keggIDtmp2,$koList);
												unless($ko{$koID_tmp}){
													$ko{$koID_tmp} = 1;
												}
												if($ecID_tmp){
													unless($ec{$ecID_tmp}){
														$ec{$ecID_tmp} = 1;
													}
												}
											}
										} else {
											unless($keggID =~ /$keggIDtmp\;/){$keggID .= $keggIDtmp.";";}
											my ($koID_tmp,$ecID_tmp) = getKO($keggIDtmp,$koList);
											unless($ko{$koID_tmp}){
												$ko{$koID_tmp} = 1;
											}
											if($ecID_tmp){
												unless($ec{$ecID_tmp}){
													$ec{$ecID_tmp} = 1;
												}
											}
										}
									} elsif($_ =~ /orf(.)+/){	### only for canal
										my $keggIDtmp = $&;
										$keggIDtmp =~ s/orf/CaO/;
#		print $keggIDtmp,"\n";
										unless($keggID =~ /$keggIDtmp\;/){
											$keggID .= $keggIDtmp.";";
#		print "HERE:$keggID\n";
											my ($koID_tmp,$ecID_tmp) = getKO($keggIDtmp,$koList);
											unless($ko{$koID_tmp}){
												$ko{$koID_tmp} = 1;
											}
											if($ecID_tmp){
												unless($ec{$ecID_tmp}){
													$ec{$ecID_tmp} = 1;
												}
											}
										}
									}
								}
		#						my $keggID = $uniHit_tmp[@uniHit_tmp-1];
		#						($koID,$ecID) = getKO($keggID,$koList);
							} else {
								print "Cannot find hit in $specID.uniprot for $specID !!\n";<>;
								print ERROR "$protID\t$$protID\t$specID.uniprot\n";
							#	unless($lackUniProt{$specID}){$lackUniProt{$specID} = $protID;}
							#	else {$lackUniProt{$specID} .= "\t".$protID;}
							}
							$koID = join("\t",keys(%ko));	$koID =~ s/\t\t+/\t/g;	$koID =~ s/^\t//; $koID =~ s/\t$//;
							$ecID = join("\t",keys(%ec));	$ecID =~ s/\t\t+/\t/g;	$ecID =~ s/^\t//; $ecID =~ s/\t$//;
#		print "KO=$koID\n";
#	}
						}

						### other species
						else {
							$keggID = $protID;
							($koID,$ecID) = getKO($keggID,$koList);
	#						if(length $koID > 2) {print "here: $koID";<>;}
						}

						### print KO (and EC if available)
						if(length($koID) > 2){
	#						print $koID," - ",$ecID;<>;
							print OUT "$specID:$protID","\t",$koID,"\t",$ecID,"\t",$keggID,"\n";
						}
#						 else {
#							print $protID," - ",$keggID," has no KO!\n";<>;
#						}
					} else {
					#	print "NO KO FILE FOR THIS SPECIES $specID\n";<>;
					}
				}
			}
			close (OUT);
		}
	}
}

exit;

sub getKO {
	my ($keggID,$koList) = @_;
	my $koID = ""; my $ecID = "";
	if($koList =~ /\n$keggID\b(.)*?\n/){
		my $hit = $&; chomp($hit);
		my @koHit = split(/\t/,$hit);
		$koID = $koHit[1];
		if($koHit[2]){
			$ecID = $koHit[2]; $ecID =~ s/\n//g;
		}
	#	print $koID," - ",$ecID,"\n";
	}
	return ($koID,$ecID);
};
