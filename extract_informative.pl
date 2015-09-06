#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;
use Bio::SeqIO;

my %opt;
my $usage =  "USAGE:
	$0 
	-s sam file of short reads aligned to reference  if not specifed , it will read from stdin 
	-n the ID of TE seq 
	-p out_put file prefix
	-h print this help
	-g genome file
	" ;


die $usage if ( @ARGV == 0);
getopts("s:n:g:p:h",\%opt);
die $usage if ($opt{h});



my $id = $opt{n};
my $pre = $opt{p};
my $sam = $opt{s};
my $genome = $opt{g};

open my $fh,">$pre.$id.informative.sam" or die $!;
my $seq_in = Bio::SeqIO -> new ( -file => $genome,-format=>"fasta");
my %gid;
while(my $seq = $seq_in -> next_seq){
	my $id = $seq -> id;
	$gid{$id} =1 ;
}
my %reads;      # hash to assistant to detect paired reads in sam file
my @rs;			# array contain the pairs of one fragment

open SAM, $sam or die $!;


my $num_infor =0;
while (<SAM>){     #  reading sam file one by one 
	chomp;
	if (/^@/){     #   save header 
		print $fh "$_\n"; 
		next;
	}
	
	my ($title,$chr,$cig,$rnext) = (split /\t/,$_)[0,2,5,6];
	
	if (! keys(%reads )  or exists ($reads{$title})){     #if have no defined hash %reads or the key of %reads is equal to the $title, then I am reading another pair read
		push @rs,$_;
		$reads{$title} = "1";
		if(eof(SAM)){
			print_clu(@rs);
		}
	}else{
		
		print_clu(@rs);

		undef(%reads);
		$reads{$title} = 1;
		@rs = ($_);
		print_clu(@rs) if (eof(SAM));
	}
}

if($num_infor == 0){
	die "NO insertions can be found. Exit!\n";
}

sub print_clu{
	my @va = @_;
	my $pt = join "\n",@va;     # $pt is ready to print
	#my $map_q;     # used to check map_q of reads at genome
	
	my $gboo;
	my $tboo;
	foreach my $it (@va){
		my ($title,$flag,$chr,$mq,$cig,$rnext) = (split /\t/,$it)[0,1,2,4,5,6];
		return unless ( $gid{$chr} or  $chr eq $id);
		return unless ( $gid{$rnext} or  $rnext eq $id or $rnext eq "=" );
		if($gid{$chr} or $gid{$rnext}){
			$gboo = 1;
		}
		if ($chr eq $id or $rnext eq $id){
			$tboo = 1;
		}
	}			
	if ($gboo and $tboo){
		print $fh "$pt\n" ;
		$num_inf ++;
	}
}

