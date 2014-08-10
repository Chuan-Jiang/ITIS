#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	my($chr,$s,$e,$d,$type,$head,$tail) = (split /\t/,$_)[0,1,2,4,5,12,13];

	next unless ($s =~ /\d/);
	my $boo = 1;
	if ($type ne "1p1"){
		$boo =0;
	}elsif($head+$tail <3){
		$boo = 0;
	}

	$s --;
	$d = ($d=~ /antisense/)?"-":"+";
	print "$chr\t$s\t$e\tTEMP\t.\t$d\n" if ($boo);
}

