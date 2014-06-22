#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Bio::SeqIO;
use Getopt::Std;

my %opt;
getopts("g:b:t:ho:l:",\%opt) ;

die "USAGE $0 
	-g genome seq file
	-t te seq file
	-h help
	-o out_put_file 
	-l name of file conating pre existes homelogus_seq
	" if ( $opt{h});

my $seq_in = Bio::SeqIO -> new (-file => $opt{g},-format => "fasta");
my %genome ;
my @order;
while (my $seq_obj = $seq_in -> next_seq){
	my $id = $seq_obj -> id;
	my $seq = $seq_obj -> seq;
	$genome{$id} = $seq;
	push @order,$id;
}
my %te = Seq::seq_hash($opt{t});

print STDERR "making index for blastn\n";
system ("makeblastdb -in=$opt{g} -dbtype='nucl'") == 0 or die $!;
open BLA, "blastn -query $opt{t} -outfmt 6 -db $opt{g}|" or die $!;

open LS, ">$opt{l}" or die $!;
while(<BLA>){
	chomp;
	print LS "$_\n";
	my ($chr,$s,$e) = (split /\t/,$_)[1,8,9];
	($s,$e) = sort {$a<=>$b}($s,$e);
	my $l = $e-$s+1;
	#print STDERR "substracting...\nalignment $_\n";
	substr($genome{$chr},$s-1,$l) = "N"x$l;
}

open OUT, ">$opt{o}" or die $!;
foreach my $k ( @order){
	print OUT ">$k\n$genome{$k}\n";
}
my ($te_n,$te_seq) = each (%te);
print OUT ">$te_n\n$te_seq\n";


