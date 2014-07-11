#!/usr/bin/perl
use warnings; use strict;
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


open my $fh, $file or die $!;
my $tmp = "tmp_bwa".time();
open my $out, ">$tmp" or die $!;

my $exi = 1;
while (<$fh>){
	chomp;
	my ($id,$r,$chr,$seq)  = (split /\t/,$_)[0,1,2,9];	
	#$r =~ s/.*(\d)/$1/;
	
	if ($chr =~ /$te/  ){
		my $q = "J"x(length($seq));
		print $out "\@$id:$exi\n$seq\n\+\n$q\n";
		$exi ++;
	}
}

my $o = $opt{p};

system ("bwa mem -T 20 $opt{i} $tmp -a >$o.alnte.sam" )== 0 or die $!;


system("rm -rf $tmp") == 0 or die $!;

