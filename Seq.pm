package Seq;
use warnings; use strict;
use Bio::SeqIO;
use Bio::Seq;

sub seq_hash{
	my $file = shift @_;
	my $seq_in = Bio::SeqIO -> new (-file => $file,-format => "fasta");
	my %hash;
	while (my $seq_obj = $seq_in -> next_seq){
		my $id = $seq_obj -> id;
		my $seq = $seq_obj -> seq;
		$hash{$id} = $seq;
	}
	return %hash;
}

sub rev_com{
	my $seq = shift @_;
	(my $seq_com = $seq )=~ tr/ATCGatcg/TAGCtacg/;
	my $rev_com = reverse $seq_com;
	return  $rev_com;
}

sub translate{
	my $seq = shift @_;
	my $seq_obj = Bio::Seq -> new (-alphabet=>"dna",-seq => $seq);	
	my $pro = $seq_obj -> translate;
	my $pro_seq = $pro -> seq;
	return $pro_seq;
}
1;

