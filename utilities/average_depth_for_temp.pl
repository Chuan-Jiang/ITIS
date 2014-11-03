#!/usr/bin/perl
use warnings; use strict;

my $dep = 0;
while(<>){
	chomp;
	my @ar = split /\t/;
	my ($chr,$s,$e,$p,$d) = @ar[0,1,2,6,7];
	my $len = $e - $s;
	if($p == $len){
		$dep += $d;
		my $ever = $dep/$p;
		print "$chr:$s-$e\t$ever\n";
		$dep = 0;
	}else{
		$dep += $d;
	}
}




