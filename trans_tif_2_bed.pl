#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	my($chr,$start,$end,$res) = split(/\t/,$_,4);
	$res =~ s/\t/,/g;
	my $d;
	if ($res =~ /forward/){
		$d = "+";
	}else{
		$d = "-";
	}

	($start,$end) = sort{$a<=>$b} ($start,$end);
	$start --;
	print "$chr\t$start\t$end\tTIF,$res\t.\t$d\n";
}

