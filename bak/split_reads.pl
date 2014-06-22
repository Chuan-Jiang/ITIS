#!/usr/bin/perl
use warnings; use strict;
use Seq;

while (<>){
	chomp;
	#open my $r1,">>reads1.fq" or die $!;
	#open my $r2,">>reads2.fq" or die $!;
	if ($_ =~ /^@/){
		next;
	}
	my ($id,$flag,$cig,$seq,$qal) = (split /\t/,$_)[0,1,5,9,10];
	next if ($cig =~ /H/);
	if ($flag =~ /r/){
		$seq =  Seq::rev_com($seq);
		$qal = reverse $qal;
	}
	
	my $len = length $seq;
	(my $r = $flag) =~ s/.*(\d)/$1/;
	#my $fh = ($r ==1?$r1:$r2);

	for my $i ( 0..($len-20)){
		my $s = substr($seq,$i,20);
		my $q = substr($qal,$i,20);
		print  "\@$id#$r#$i\n$s\n+\n$q\n";
	}
}


