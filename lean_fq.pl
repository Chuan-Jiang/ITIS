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
my $bwa = "bwa";
open R1,">$pre.fq1" or die $!;
open R2, ">$pre.fq2" or die $!;

open ALN, "$bwa mem  -MT 20 -t $cpu $index $r1 $r2 2>/dev/null | samtools view -XS - |" or die $!;
my %reads;
my $last ;

while (<ALN>){
	chomp;
	next if (/^@/);

	my($id,$flag,$seq,$q) = (split /\t/)[0,1,9,10];
	
	if($last and $id ne $last){
		prt($last);
	}
	if(eof(ALN)){
		prt($id);
	}
	

	next if ($flag =~ /uU/);
	next if ($flag =~ /s/);
	my $r = $flag =~ /1/? 1:2;

	if($flag =~ /r/){
		$seq = reverse($seq);
		$seq =~ tr/atcgATCG/tagcTAGC/;
	}
	$reads{$id}{$r} =">$id\n$seq\n";
	$last = $id;
}
	
sub prt{
	my $id = shift @_;
	my $r1 = $reads{$id}{1};
	my $r2 = $reads{$id}{2};
	if($r1 and $r2){
		print R1 $r1;
		print R2 $r2;
	}
	undef $reads{$id};
}





