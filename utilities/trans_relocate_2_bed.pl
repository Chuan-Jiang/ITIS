#!/usr/bin/perl
use warnings; use strict;


while (<>){
	chomp;
	next if (/^#/);

	my( $chr,$s,$e,$d) = (split /\t/,$_)[0,3,4,6];
	$s--;
	print "$chr\t$s\t$e\tRelocaTE\t.\t$d\n";

}
