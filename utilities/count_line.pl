#!/usr/bin/perl
use warnings; use strict;

my @files = @ARGV;

foreach my $f(@files){
	open my $hd, $f or die $!;

	my @lines = <$hd>;
	print "$f\t".scalar @lines."\n";
}

