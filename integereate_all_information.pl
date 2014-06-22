#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Getopt::Std;

my %opt;
getopts("b:a:g:i:t:h",\%opt);

my($bed_file,$gff_file,$genome_file,$tnt) = ($opt{b},$opt{a},$opt{g},$opt{t});
die  "USAGE:
	$0 :
	-b bed file of insertion locs
	-a annotation file in gff format ( only contain 'GENE')
	-g genome sequence
	-i generate igv snpshot batch file? <file name> 
	-t the name of your retrotransposon name
	-h print this help infor
	" if ($opt{h});

################ generate igv batch   #############
if ($opt{i}){
	open FH,">$opt{i} " or die $!;
	open BED,"$bed_file" or die $!;
	print FH "snapshotDirectory ./$opt{t}_snap\n";
	while (<BED>){
		chomp;
		 my ($chr,$s,$e) = split /\t/,$_;
		 if (($e - $s) > 300){
			 die "$chr,$s,$e,haha\n";                                                                                        
		 }else{
			 my $slop = 150 - ($e - $s);
			 my $s_p = $s - int($slop/2);
			 my $e_p = $e + int($slop/2);
			 print FH "goto $chr:$s_p-$e_p\n";
			 print FH "snapshot ${chr}_${s}_${e}_slop$slop.png\n";
		 }
	 }
 }



################# generate usefull gff file containg intergenic region  ##########
open TEM,">tem" or die $!;
open GFF,"$gff_file " or die $!;
chomp (my @gff = <GFF>);
for (my $i = 1;$i<@gff;$i++){
	print TEM "$gff[$i-1]\n";
	my $l = $gff[$i-1];
	my $n  = $gff[$i];
	my @las = split /\t/,$l;
	my @nos = split /\t/,$n;
	if (($nos[3] - $las[4]) > 0 and $nos[0] eq $las[0]){
		my $s = $las[4] +1;
		my $e = $nos[3] -1;
	
		my %la;
		foreach  (split /;/,$las[8]){
			my ($k,$v) = split /=/,$_;
			$la{$k} = $v;
		}

		my %ne;
		foreach (split /;/,$nos[8]){
			my ($k,$v) = split /=/,$_;
			$ne{$k} = $v;
		}

		my $la_ge = $la{ID};
		my $la_no = $la{Note};

		my $ne_ge = $ne{ID};
		my $ne_no = $ne{Note};
		print TEM "$las[0]\t$las[1]\tIntergenic\t$s\t$e\t.\t.\t.\tID=${la_ge}_${ne_ge};Note=Intergenic_${la_no}_INSERT_$ne_no\n";
	}
}
##################### put genome in hash##################################

my %genome = Seq::seq_hash($genome_file);
##########################################################################



################### intersect using bedtools  ###########################

open BEDTOOL, "bedtools intersect -a $bed_file -b tem -wa -wb |" or die $!;

my $la=0;
my $la_t;
while (<BEDTOOL>){
	chomp;
	#print "$_\n";
	my($chr,$s,$e,$d,$t,$an) = (split /\t/,$_)[0,1,2,3,6,12];
	next if (($s+$e) == $la and $la_t eq "gene");
	my ($up,$down);
	my $tnt_seq = ($d =~ /:S:/)?$genome{$tnt}:Seq::rev_com($genome{$tnt});
	if ($d =~ /:P$/){
		$up = substr ($genome{$chr},$e-1005,1005);
		$down = substr ($genome{$chr},$s,1005);
	}else{
		$up = substr ($genome{$chr},$s-1000,1000)."=";
		$down = "=".substr ($genome{$chr},$e-1,1000);
	}
	my $seq = $up.lc($tnt_seq).$down;
	print "$chr\t$s\t$e\t$d\t$t\t$an\t$seq\n";
	$la = $s+$e;
	$la_t = $t;

}

unlink "tem";
