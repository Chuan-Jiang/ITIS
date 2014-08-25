#!/usr/bin/perl
use warnings; use strict;
use Bio::SeqIO;
# this scripts use to randomly insert 10000; to calculate the distributuon
#
#
my $genome_file = shift @ARGV;
my $seq_in = Bio::SeqIO -> new (-file => $genome_file, -format => "fasta");

my @chrs;
my %geno;
my $tot;
while(my $seq_obj = $seq_in -> next_seq()){
	my $id = $seq_obj -> id();
	my $len = $seq_obj -> length();
	
	push @chrs,$id;
	$tot += $len;
	$geno{$id} = $tot;
}


srand(100);
my @nums;
for (1..10000){
	push @nums ,int(rand($tot-5));
}
@nums = sort {$a <=> $b} @nums;



my $pp = 0;

LOOP: foreach my $id (@chrs){

	while(1){
		if(@nums == 0){
			last;
		}elsif($nums[0] <=  $geno{$id}-5){
			my $first = $nums[0] - $pp;
			my $last = $first + 5;
			print "$id\t$first\t$last\n";
			shift @nums;
		}else{
			$pp = $geno{$id}- 1;
			next LOOP;
		}
	}
}





