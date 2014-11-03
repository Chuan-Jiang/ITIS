#!/usr/bin/perl
use warnings; use strict;

my $file = shift @ARGV;


open GENE, "bedtools window -w 0 -a ~/jiangchuan/genome/medicago/Mt4.0v1_genes_20130731_1800.gff3 -b $file | " or die $!;

my %hits;

while (<GENE>){
	my@ar = split /\t/;
	if ($ar[2] =~ /exon/){
		$hits{$ar[9],$ar[10],$ar[11]}{0} = "exon";
		next;
	}elsif($ar[2] =~ /gene/){
		$hits{$ar[9],$ar[10],$ar[11]}{0} = "gene";
	}
}

count_hit("-r 2500 -l 0","r25");
count_hit("-l 2500 -r 0 ","l25");
#count_hit("-l 5000","l50");
#count_hit("-r 5000","r50");


my ($exon,$intron,$l25,$r25,$inter) ;
open RE, "$file" or die $!;
while(<RE>){
	my@ar = split /\t/;
	my $va = $hits{$ar[0],$ar[1],$ar[2]};
	unless($va){
		$inter ++;
	}else{
		my %h = %$va;
		#print keys %h;
		#print ": $ar[0],$ar[1],$ar[2]\n";
		foreach my $k ( sort {$a <=> $b} keys %h){
			my $v = $h{$k};
			if($v =~ /exon/){
				$exon ++;
			}elsif($v =~ /gene/){
				$intron ++;
			}elsif($v =~ /l25/){
				$l25 ++;
			}elsif($v =~ /r25/){
				$r25 ++;
			}
			last;
		}
	}
}
print "$exon,$intron,$l25,$r25,$inter\n";


sub count_hit{
	my($pa,$v) = @_;
	open IN , "bedtools window $pa -a ~/jiangchuan/genome/medicago/Mt4.0v1_genes_20130731_1800.gff3 -b $file |" or die $!;
	while (<IN>){
		my $dis;
		my@ar = split /\t/;
		if($ar[2] =~ /gene/ and ! exists $hits{$ar[9],$ar[10],$ar[11]}{0} ){
			my $h = $ar[3] - $ar[11];
			my $t = $ar[10] - $ar[4];
			$dis = $h>0?$h:$t;
			die if ($dis <0);
			$hits{$ar[9],$ar[10],$ar[11]}{$dis}= "$v";
		}
	}
}
