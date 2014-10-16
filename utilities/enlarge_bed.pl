#!/usr/bin/env perl
use warnings; use strict;

die  "$0 bed_file window_size\n" unless @ARGV;

open BED, shift @ARGV or die $!;
my $w = shift @ARGV;
while(<BED>){
	chomp;
	my($chr,$s,$e,$rest) = split (/\t/,$_,4);
	my $dis = $e-$s+1;
	if($dis<$w){
		my $pad = int(($w-$dis+1)/2);
		$s -= $pad;
		$e += $pad;
	}
	print "$chr\t$s\t$e\t$rest\n";
}
print "chr1\t1\t4\t.\t.\t.\n";

