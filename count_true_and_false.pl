#!/usr/bin/perl
use warnings; use strict;

my($ref_list,$list) = @ARGV;

my $ref_ln;
open REF, $ref_list or die $!;
$ref_ln ++ while(<REF>);

my $ln;
open LIS, $list or die $!;
$ln++ while (<LIS>);


my($t_500,$e_500) = count(500);
my($t_100,$e_100) = count(100);

print "$ref_ln\t$ln\t$t_500\t$t_100\t$e_100\n";

sub count{
my $w  = shift @_;

open WIN, " bedtools window -a $ref_list -w $w -b $list | " or die $!;

my $tot_over = 0;
my $tot_exact = 0;
my $rec=0;
while (<WIN>){
	 chomp;
	 my($s_r,$e_r,$d_r,$s,$e,$d) = (split /\t/,$_)[1,2,3,5,6,9];
	 if($s_r != $rec){
		 my $sign = ($d_r == 1)? "\\+":"\\-";
		 if ($d =~ /^$sign/ or $d =~ /NA|\./){
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
return("$tot_over","$tot_exact");
}





