#!/usr/bin/perl
use warnings; use strict;

while (<>){
	my ($chr,$s,$e,$dir) = split /\t/;
	if ($e - $s == 5){
	my 	( $d ) = $dir =~ /^(\+|-)/;
		$s -= 200;
		$e += 200;
		print "$chr\t$s\t$e\t.\t.\t$d\n";
	}
}

