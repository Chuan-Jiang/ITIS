#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Getopt::Std;

my %opt;

my $usage =  "USAGE $0
	-n te seq ID in fasta
	-s full informative sam of TE
	-i full path bwa index of te sequence
	-h print this help
	-p prefix of you output file
	" ;

die "$usage\n"  if (@ARGV == 0);
getopts("n:s:i:p:h",\%opt);
die "$usage\n" if ($opt{h});

my ($te,$file) = ($opt{n},$opt{s});
my $bwa = "bwa";

open my $fh, $file or die $!;
my $tmp = "tmp_bwa".time();
open my $out, ">$tmp" or die $!;

while (<$fh>){  #  open full informative sam file
	chomp;
	next if (/^@/);
	my ($id,$flag,$chr,$seq)  = (split /\t/,$_)[0,1,2,9];	
	(my$r) = $flag  =~ /(\d)/;
	print "ERRO:::$id\n" unless $seq;
	$seq = Seq::rev_com($seq) if ( $flag =~ /r/);
	#print "$id\t$seq\n" if  ( $flag =~ /r/);
	if ($chr =~ /$te/  ){
		my $q = "J"x(length($seq));
		print $out "\@$id:$r\n$seq\n\+\n$q\n";
	}
}

my $o = $opt{p};

system ("$bwa mem  -T 20 $opt{i} $tmp -a 2>/dev/null >$o.$te.alnte.sam" )== 0 or die $!;

system("rm -rf $tmp") == 0 or die $!;

