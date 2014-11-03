#!/usr/bin/perl
use warnings; use strict;

my($ins_bed,$sup_bam)= @ARGV;

open INS, $ins_bed or die $!;
while(<INS>){
	chomp;
	my ($chr,$s,$e,$ann) = split /\t/;
	$e++;

	if($ann =~ /GT=NA/){
		$s = $s - 100 ;
		$e = $e + 100;
	}

	### find the forward primer
	open SAM, "samtools view $sup_bam $chr:$s-$e |" or die $!;
	print "samtools view $sup_bam $chr:$s-$e \n";	
	my @reads = <SAM>;
	unless (@reads){
		print "$_\tNA\tNA\n";
	}else{
		my $fir = shift @reads;
		my $las = pop @reads;
		
		my ($l5,$p5) = ext_seq(5,$fir);
		my ($l3,$p3) = ext_seq(3,$las);
		
		if($l5 >= $s){
			$p5 = "NA" ;
		}
		
		if($l3 <= $s){
			$p3 = "NA";
		}else{
			$p3 = reverse($p3);
			$p3 =~ tr/ATGC/TACG/;
		}
		print "$_\t$p5\t$p3\n";
	}
}

sub ext_seq{
	my ($ty,$aln) =  @_;
	my @arr = split /\t/, $aln;
	my $pos = $arr[3];
	my $cig = $arr[5];
	my $seq = $arr[9];
	my $seq_r = $seq;
	if ($ty == 5 and $cig =~ /(\d+)S$/){
		$seq_r = substr($seq,0,-$1);
	}elsif($ty == 3 and $cig =~ /^(\d+)S/){
		$seq_r = substr($seq,$1-1);
	}
	return ($pos,$seq_r);
}
