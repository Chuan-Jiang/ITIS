#!/usr/bin/perl
use warnings; use strict;

unless (-e "snap_shot"){
	print "snapshotDirectory ./snap_shot\n" ;
}else{
	die "dir exists\n";
}

while (<>){
	chomp;
	my ($chr,$s,$e) = split /\t/,$_;
	if (($e - $s) > 300){
		die "$chr,$s,$e,haha\n";
	}else{
		my $slop = 150 - ($e - $s);
		my $s_p = $s - int($slop/2);
		my $e_p = $e + int($slop/2);
		print "goto $chr:$s_p-$e_p\n";
		print "snapshot ${chr}_${s}_${e}_slop$slop.png\n";
	}
}
