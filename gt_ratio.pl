#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	if($_ =~ /GT=(\d+),(\d+)/){
		print "$1\t$2\n";
	}
}

