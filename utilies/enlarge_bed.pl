#!/usr/bin/env perl
use warnings; use strict;


open BED, shift @ARGV or die $!;
my $w = shift @ARGV;
while(<BED>){
	chomp;
	my($chr,$s,$e,$rest) = split (/\t/,$_,4);
	$s -= $w;
	$e += $w;
	print "$chr\t$s\t$e\t$rest\n";
}
print "chr1\t1\t4\t.\t.\t.\n";

