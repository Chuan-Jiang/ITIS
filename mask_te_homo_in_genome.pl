#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Bio::SeqIO;
use Getopt::Std;

my %opt;
getopts("g:t:o:h",\%opt) ;

die "USAGE $0 
	-g genome seq file
	-t te seq file
	-o out_put_file 
	-h help
	" if ( $opt{h});

# put genome seq in hash
my $seq_in = Bio::SeqIO -> new (-file => $opt{g},-format => "fasta");
my %genome ;
my @order;
while (my $seq_obj = $seq_in -> next_seq){
	my $id = $seq_obj -> id;
	my $seq = $seq_obj -> seq;
	$genome{$id} = $seq;
	push @order,$id;
}
# put te seq in hash
my %te = Seq::seq_hash($opt{t});

# using blast2seq to identify the te homolog
open BLA, "blastn -query $opt{t} -subject $opt{g} -outfmt 6 |" or die $!;
while(<BLA>){
	chomp;
	my ($chr,$s,$e) = (split /\t/,$_)[1,8,9];
	($s,$e) = sort {$a<=>$b}($s,$e);
	my $l = $e-$s+1;
	#print STDERR "substracting...\nalignment $_\n";
	substr($genome{$chr},$s-1,$l) = "N"x$l;
}



# put masked genome seq and te seq together 
open OUT, ">$opt{o}" or die $!;
foreach my $k ( @order){
	print OUT ">$k\n$genome{$k}\n";
}

foreach my $k(keys %te){
	print OUT ">$k\n$te{$k}\n";
}
