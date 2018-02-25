#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
check connectivity of annotated yeast proteins
30.11.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl checkConnect.pl -i sacce_unk.list.NEW.KO -p yeast.ppi -k /share/project/vinh/pathways_kegg/keggDB -m KEGG_annotation/refseq/sacce.refseq\n";
    print "-i\tinput KO annotated file\n";
    print "-p\tyeast ppi file\n";
    print "-k\tKEGG database folder\n";
    print "-m\tyeast refseq file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_p,$opt_k,$opt_m);
getopts('i:p:k:m:');

my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input annotated file given\n");
my $ppiFile = ($opt_p) ? $opt_p : usage("ERROR: No PPI file given\n");
my $kegFol = ($opt_k) ? $opt_k : usage("ERROR: No KEGG database folder given\n");
my $refseqFile = ($opt_m) ? $opt_m : usage("ERROR: No refseq file given\n");

### read and parse KEGG file
open(KO,"$kegFol/ko_map.list") || die "Cannot open $kegFol/ko_map.list!\n";
my %ko2path; # $path{koID} = all_pathIDs
foreach my $line(<KO>){
  chomp($line);
  my @tmp = split(/\t/,$line);  # K00001	ko00010
  unless($ko2path{$tmp[0]}){
    $ko2path{$tmp[0]} = $tmp[1];
  } else {
    $ko2path{$tmp[0]} .= ";".$tmp[1];
  }
}
close (KO);

open(MAP,"$kegFol/map_desc.list") || die "Cannot open $kegFol/map_desc.list!\n";
my %pathType; # $pathType{$pathID} = path_type
foreach my $line(<MAP>){
  chomp($line);
  my @tmp = split(/\t/,$line);  # 00020	Citrate cycle (TCA cycle)	Carbohydrate metabolism	Metabolism
  $pathType{"ko$tmp[0]"} = $tmp[2];
}
close (MAP);

### read and parse PPI file
open(PPI,$ppiFile) || die "Cannot open $ppiFile!\n";
my %ppi;  # $ppi{seqID} = all_connected_ids
foreach my $line(<PPI>){
  chomp ($line);
  my @tmp = split(/\t/,$line);
  unless($ppi{$tmp[0]}){
    $ppi{$tmp[0]} = $tmp[1];
  } else {
    $ppi{$tmp[0]} .= ";".$tmp[1];
  }
}
close (PPI);

### read and parse refseq file to map refseq IDs with in-house IDs
open(REF,$refseqFile) || die "Cannot open $refseqFile\n";
my %mapID;  # $mapID{inhouseID} = refseqID
foreach my $line(<REF>){
  chomp($line);
  my @tmp = split(/\t/,$line);
  my $refID = $tmp[0]; $refID =~ s/sacce_2336_1://;
  my $inhouseID = "sacce:".$tmp[1];
  $mapID{$inhouseID} = $refID;
}
close (REF);

### check connectivity for annotated proteins
open(IN,$inFile) || die "Cannot open $inFile\n";
my @in = <IN>;
close (IN);
my $in = join("",@in);
open(PPI,">$inFile.PPI");
print PPI "protID\trefseqID\tnodeDegree\n";
open(KO2PATH,">$inFile.ko2path");
print KO2PATH "KO\tnb_pathways\tpathwayID\tpathwayType\n";

foreach my $line(@in){
  chomp($line);
  if($line =~ /K\d{5}/){
    my @tmp = split(/\t/,$line);  # sacce:2037	K08870;0.99956386422;dme:Dmel_CG3172	K14962;0.99634282393;dre:321614
    my $seqID = shift(@tmp);

    # check ppi
    print PPI $seqID,"\t",$mapID{$seqID};
    if($ppi{$mapID{$seqID}}){
      my @partners = split(/;/,$ppi{$mapID{$seqID}});
      # print OUT "\t",scalar(@partners),"\t",$ppi{$mapID{$seqID}};
      print PPI "\t",scalar(@partners);
    } else {
      # print OUT "\t0\tNA";
      print PPI "\t0";
    }
    print PPI "\n";

    # check pathway
    foreach my $hit(@tmp){
      my $type = "";
      my @hitTMP = split(/;/,$hit); # K08870;0.99956386422;dme:Dmel_CG3172
      if($ko2path{$hitTMP[0]}){
        my @pathways = split(/;/,$ko2path{$hitTMP[0]});
        foreach my $path(@pathways){
          unless($type =~ /$pathType{$path}/){
            $type .= ";".$pathType{$path};
          }
        }
        $type =~ s/^;//;
        print KO2PATH $hitTMP[0],"\t",scalar(@pathways),"\t",$ko2path{$hitTMP[0]},"\t",$type,"\n";
      } else {
        print KO2PATH "$hitTMP[0]\t0\tNA\tNA\n";
      }
    }
  }
}
close (PPI);
close (KO2PATH);
print "FINISHED! Check output $inFile.PPI and $inFile.ko2path\n";
exit;
