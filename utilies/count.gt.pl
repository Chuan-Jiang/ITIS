#!/usr/bin/perl
use warnings; use strict;

my ($su,$nsu);
while(<>){
	chomp;
	my($aa,$bb) = $_ =~ /GT=(\d+),(\d+)/;
	$su += $aa;
	$nsu += $bb;
}
print "$su\t$nsu\n";

