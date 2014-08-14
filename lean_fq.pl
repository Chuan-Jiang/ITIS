#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

my $usage = "usage: $0
	-1 : read1
	-2 : reads
	-p : prefix of out put reads  
	-i : reference index file
	-c : cpu number for bwa
	-h : print this help

	BWA should in your PATH
	";
die $usage if (@ARGV ==0);
my %opt;
getopts("1:2:p:i:c:h",\%opt);

my($r1,$r2,$pre,$index,$cpu) = ($opt{1},$opt{2},$opt{p},$opt{i},$opt{c});

open R1,">$pre.fq1" or die $!;
open R2, ">$pre.fq2" or die $!;

open ALN, "bwa mem  -T 20 -t $cpu $index $r1 $r2 | " or die $!;
my %reads;
my $last ;

while (<ALN>){
	chomp;
	next if (/^@/);

	my($id,$flag,$seq,$q) = (split /\t/)[0,1,9,10];
	
	if($last and $id ne $last  and %reads ){
		print  R1 $reads{1};
		print  R2 $reads{2};
		undef(%reads);
	}
	$last = $id;

	$flag = unpack("B32",pack("N",$flag));
	my @ar = (split //,reverse($flag))[0..11];
	next if ($ar[2] and $ar[3]);
	next if ($ar[11] or $ar[8]);
	my $r = $ar[6]?1:2;
	die "lean fq error\n" if ($r == 2 and $ar[7] == 0);
	if($ar[4] == 1){
		$seq = reverse($seq);
		$seq =~ tr/ATCG/TAGC/;
	}
	$reads{$r} ="\@$id\n$seq\n+\n$q\n";
}
	
	
	


