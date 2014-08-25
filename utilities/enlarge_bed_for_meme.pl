#!/usr/bin/perl
use warnings; use strict;
my($file,$step) = @ARGV;

open my $fh, "$file " or die $!;

while (<$fh>){
	my ($chr,$s,$e,$dir) = split /\t/ , $_, 4 ;
	if ($e - $s == 5){
		$s -= $step;
		$e += $step;
		print "$chr\t$s\t$e\t$dir";
	}
}

