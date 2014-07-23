#!/usr/bin/perl
use warnings; use strict;

my($ref_list,$list) = @ARGV;

my $ref_ln;
open REF, $ref_list or die $!;
$ref_ln ++ while(<REF>);

my $ln;
open LIS, $list or die $!;
$ln++ while (<LIS>);




open WIN, " bedtools window -a $ref_list -w 100 -b $list | " or die $!;



my $tot_over = 0;
my $tot_exact = 0;
my $rec=0;
while (<WIN>){
	 chomp;
	 my($s_r,$e_r,$d_r,$s,$e,$d) = (split /\t/,$_)[1,2,3,5,6,7];
	 if($s_r != $rec){
		 my $sign = ($d_r == 1)? "\\+":"\\-";
		 if ($d =~ /^$sign/ or $d =~ /NA/){
			 $tot_over ++;
			 if($s_r == $s and $e_r == $e){
				 $tot_exact ++;
			 }
		 }
	 }else{
		 next;
	 }
	 $rec = $s_r;
 }

 print "$ref_ln\t$ln\t$tot_over\t$tot_exact\n";



