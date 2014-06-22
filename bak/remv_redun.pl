#!/usr/bin/perl
use warnings; use strict;

my($ids,$fq)= @ARGV;

open ID,$ids or die $!;
my @ids ;
while(<ID>){
	chomp;
	$_ =~ s/\s*(\S+)\s*/$1/;
	push @ids,$_;

}

open FQ,$fq or die $!;
while (<FQ>){
	chomp;
	foreach my $i(@ids){
		if($_ =~ /$i/){
			<FQ>;
			<FQ>;
			<FQ>;
		}
	}
	print "$_\n";
}
