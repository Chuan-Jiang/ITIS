#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

my $help = "$0 
	-i : insertion list in bed format
	-l : postion list of te homo in ref genome
	-n : default <2,3,1,1>, the minimum requried reads supporting insertion, minimum required reads cover TE start siteand and minimum required reads cover TE end site
	-q : degault <1>, the minimum required average mapping quality
	-d : default <2-200>, the reads depth range
	-h : help message
	";

die $help unless ( @ARGV);

my %opt;
getopts("i:l:c:n:q:d:r:h",\%opt);
die $help  if($opt{h});

######## parameters ###########

my $ins_file = $opt{i};
my $lst = $opt{l};

my ($CF,$total_reads,$te_s,$te_e) = $opt{n}?(split ',',$opt{n}):(split ',','2,3,1,1');
my $MQ = $opt{q}?$opt{q}:1;
my ($DP_L,$DP_H) = $opt{d}?(split ',',$opt{d}):(split ',',"2,200");

open INS, "$ins_file" or die $!;


## all homo in hash %homos
my %homos;
if($lst){
	open LST, $lst or die $!;
	while(<LST>){
		chomp;
		my($chr,$s,$e,$te) = split /\t/;
		$s = $s - 100;
		$e = $e + 100;
		foreach my $i ($s..$e){
			$homos{$te}{$chr}{$i} = 1;
		}
	}
}


while(<INS>){
	
	my $boo = 1;

	my($chr,$s,$e,$t,$rest) = (split /\t/, $_,5);
	
	my@tags = split /;/, $t;



### parse the rest key and values	
	my %other;
	foreach my $r (@tags){
		my($k,$v) = split /=/,$r;
		$other{$k} = $v;
	}

# filter support reads number
	if(exists $other{SR}){
		my($cf,$tot,$r1,$r2,$r3,$r4) = split /,/,$other{SR};
		if($cf < $CF or $tot < $total_reads or $r1+$r3 < $te_s or $r2+$r4 < $te_e){
			$boo = 0;
		}
	}
# filter eveage mapping valeu
	if(exists $other{MQ} and $other{MQ} < $MQ){
		$boo = 0;
	}

# filter bg depth
	if(exists $other{DP} ){
		my $dp = $other{DP};
		if ($dp < $DP_L or $dp > $DP_H){
			$boo = 0;
		}
	}
# mark ins near known site
	if ($lst){
		my $near = 0 ;
		for my $i ($s..$e){
			if ($i ~~ %{$homos{$other{NM}}{$chr}}){
				$near = 1;
			}
		}
		if($near){
			$t .= ";NB=Y";
		}else{
			$t .= ";NB=N";
		}
	}
	print join "\t", ($chr,$s,$e,$t,$rest) if $boo;
}










