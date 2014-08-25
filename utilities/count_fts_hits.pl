#!/usr/bin/env perl
use warnings; use strict;
use File::Basename;

my($fts_f,$bed_list_f) = @ARGV;


open BED, "$bed_list_f" or die $!;
chomp(my @bed_list = <BED>);
foreach my $bed (@bed_list){
	my $fn = basename($bed);
	system ("enlarge_bed.pl $bed >$fn.f1000.bed" ) == 0 or die $!;
	system ("bedtools getfasta -fi ~/jiangchuan/genome/medicago/JCVI.clean.fa -bed $fn.f1000.bed -fo $fn.f1000.fa") == 0 or die $!;
	open BLA, ("blastn -evalue 1e-10 -query $fts_f -subject $fn.f1000.fa -outfmt 6 | " ) or die $1;
	my %list;
	while (<BLA>){
		chomp;
		my($id) = (split /\t/,$_)[0];
		$list{$id} = 1;
	}
	my @k = keys %list;
	print $fn."\t".scalar@k."\n";
}

