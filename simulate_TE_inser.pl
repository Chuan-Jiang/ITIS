#!/usr/bin/env perl
use warnings; use strict;
use Bio::SeqIO;
use Seq;
use Bio::Seq;

my($genome_file,$te_file) = @ARGV;
my $seq_in = Bio::SeqIO -> new(-file=>$genome_file,-format=>"fasta");
my $seq_out = Bio::SeqIO -> new(-file=> ">$genome_file.sim.fa",-format=> "fasta");

my %te_ha = Seq::seq_hash($te_file);

my $seed = 10;
my $step= 8e6;
srand(100);
while(my $inseq = $seq_in -> next_seq){
	my $len = $inseq -> length;
	my $seq = $inseq -> seq;
	my $id = $inseq -> id;
	my $num = int($len/$step+1);
	
	my @ins_site;
	for(my $i=0; $i<$num;$i++){
		my $r = int(rand($len));
		my $tsd = substr($seq,$r-300,600);
		if($tsd =~ /n/i){
			redo;
		}else{
			push @ins_site,$r;
		}
	}
	
	my $pre = 0;
	my @frags;
	my $seq_p;
	foreach my $s ( sort {$a<=>$b} @ins_site){
		my $sr = (rand(1)<0.5)?"-1":"1";
		print "$id\t".$s."\t".($s+5)."\t$sr\n";
		
		my $te = $te_ha{tnt1};
		if($sr eq "-1"){
			$te = reverse($te);
			$te =~ tr/ATCG/TAGC/;
		}

		my $l = $s - $pre;
		my $f = substr($seq,0,$l+5);
		$seq = substr($seq,$l);
		$seq_p .= $f.lc($te);
		$pre = $s;
	}

	$seq_p .= $seq;	
	$inseq ->seq($seq_p);
	$seq_out -> write_seq($inseq);

}



