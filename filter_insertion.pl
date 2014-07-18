#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

my $help = "$0 
	-i : insertion list in bed format
	-n : default <3,1,1>, the minimum requried reads supporting insertion, minimum required reads cover TE start siteand and minimum required reads cover TE end site
	-q : degault <1>, the minimum required average mapping quality
	-d : default <3-200>, the reads depth range
	-r : default <0.2>, the minimum ratio of number of supporting reads over bg depth
	-h : help message
	";

die $help unless ( @ARGV);

my %opt;
getopts("i:n:q:d:r:h",\%opt);
die $help  if($opt{h});

my $ins_file = $opt{i};
my ($total_reads,$te_s,$te_e) = $opt{n}?(split ',',$opt{n}):(split ',','3,1,1');
my $map_q = $opt{q}?$opt{q}:1;
my ($min_d,$max_d) = $opt{d}?(split ',',$opt{d}):(split ',',"3,200");
my $ratio = $opt{r}?$opt{r}:0.2;

open INS, "$ins_file" or die $!;

while(<INS>){
	
	my $boo = 1;

	chomp;
	my($t) = (split /\t/, $_)[3];
	my@tags = split /;/, $t;
	
	my($tot,$num) = split /=/, $tags[0];
	if (abs($tot) < $total_reads){   ### filter 1
		$boo = 0;
	}
	shift @tags;

	my ($r1,$r2,$r3,$r4) = split /,/,$num;
	if ( $r1+$r3 < $te_s or $r2+$r4 < $te_e){
		$boo = 0;
	}

	
	my %other;
	foreach my $r (@tags){
		my($k,$v) = split /=/,$r;
		$other{$k} = $v;
	}
	
	if(exists $other{MQ} and $other{MQ} < $map_q){
		$boo = 0;
	}

	if(exists $other{DR} and eval($other{DR}) < $ratio){
		$boo = 0;
	}
	print "$_\n" if $boo;
}










