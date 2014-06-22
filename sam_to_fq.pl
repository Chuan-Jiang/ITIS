#!/usr/bin/perl
use warnings; use strict;
while(<>){
	my ($id,$seq)  = (split /\t/,$_)[0,9];
	my $q = "J"x(length($seq));
	print "\@$id\n$seq\n\+\n$q\n";
}
