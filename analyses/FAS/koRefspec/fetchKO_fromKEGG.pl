#!/usr/bin/env perl

use WWW::Mechanize;
use strict;
use Cwd;
use Getopt::Std;

### fetch KO information
### 02.09.2016

sub usage {
    my $msg = shift;
    print "example: perl fetchKO_fromKEEG.pl -i T00189 -o canal -n 1430\n";
    print "-i\tgenome ID\n";
    print "-o\tSpeices name\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_o,$opt_n);
getopts('i:o:n:');

my $genomeID = ($opt_i) ? $opt_i : usage("ERROR: No genome ID given!\n");
my $specID = ($opt_o) ? $opt_o : usage("ERROR: No output species name given!\n");
my $n = ($opt_n) ? $opt_n : usage("ERROR: No n given!\n");
## read URL
my $m = WWW::Mechanize->new();

# genes/proteins in a pathway
my $url = "http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+genome:".$genomeID;
$m->get($url);

## fetch infomation
my $page = $m->content('base_href' => undef);
my @lines = split(/\n/,$page);

open(OUT,">$specID.ko3");
my $c = 0;
foreach my $line(@lines){
	if($line =~ /www_bget?(.)+?"/){
if($c >= $n){
		### get protein ID
		my $protID = $&;
		$protID =~ s/\"//; $protID =~ s/www_bget\?//;

		### get KO from KO page
		my $KOurl = "http://www.genome.jp/dbget-bin/www_bget?".$protID;	# cal:CaO19.505
		my $n = WWW::Mechanize->new();
		$n->get($KOurl);
		my $ko = getKO($n->content('base_href' => undef));
		if(length $ko > 0){
			my @protTMP = split(/\:/,$protID);
#			print $protTMP[1],"\t",$ko,"\n";<>;
			print OUT $protTMP[1],"\t",$ko,"\n";
		}
}
	}
	$c ++; print $c," / ",scalar(@lines),"\n";
}
close (OUT);

print "FINISHED!! Check output at $specID.ko\n";
exit;

sub getDef{
	my $page = $_[0];
	my $desc = "";
	if($page =~ /<nobr>Definition<\/nobr><\/th>\n(.)+?\n/){
		my $tmp = $&;
		my @tmp = split(/hidden">/,$tmp);
		$desc = $tmp[@tmp-1];
		$desc =~ s/<br>//g; $desc =~ s/<\/div>//;
	}
	chomp($desc);
	return $desc;
}

sub getName{
	my $page = $_[0];
	my $name = "";
	if($page =~ /<nobr>Gene name<\/nobr><\/th>\n(.)+?\n/){
		my $tmp = $&;
		my @tmp = split(/hidden">/,$tmp);
		$name = $tmp[@tmp-1];
		$name =~ s/<br>//g; $name =~ s/<\/div>//;
	}
	chomp($name);
	return $name;
}

sub getKO{
	my $page = $_[0];
	my $ko = "";
	if($page =~ /KO<\/nobr><\/th>\n(.)+?\n/){
		my $tmp = $&;
		if($tmp =~ /dbget-bin\/www_bget\?ko\:K\d{5}/){
			$ko = $&;
			$ko =~ s/dbget-bin\/www_bget\?ko\://;
		}
	}
	chomp($ko);
	return $ko;
}
