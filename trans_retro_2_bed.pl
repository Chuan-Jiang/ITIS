#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	my($chr,$s,$e,$d) = (split /\t/)[0,1,2,5];
	print "$chr\t$s\t$e\t$d\n";
}

