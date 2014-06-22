#!/usr/bin/perl
use warnings; use strict;

while (<>){
	chomp;
	if (/^\@ACB/){
		print "$_\n";
	}elsif(/^\+/){
		print "$_\n";
	}else{
		my $p = substr ($_,0,20);
		print "$p\n";
	}
}
	
