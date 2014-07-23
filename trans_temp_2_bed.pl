#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	my($chr,$s,$e,$d) = (split /\t/,$_)[0,1,2,4];
	next unless ($s =~ /\d/);
	$s --;
	$d = ($d=~ /antisense/)?"-1":"+1";
	print "$chr\t$s\t$e\t$d\n";
}

