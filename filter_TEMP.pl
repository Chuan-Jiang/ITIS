#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	my($p,$s,$e) = (split /\t/, $_)[5,12,13];
	my $boo = 1;
	if ($p ne "1p1"){
		$boo = 0;
	}elsif($s+$e <3){
		$boo = 0;
	}
	
	print "$_\n" if ($boo);
}


