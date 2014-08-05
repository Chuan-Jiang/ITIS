#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

my %opt;
my $usage =  "USAGE:
	$0 
	-s sam file of short reads aligned to reference  if not specifed , it will read from stdin 
	-n the ID of TE seq 
	-p out_put file prefix
	-h print this help
	" ;


die $usage if ( @ARGV == 0);
getopts("s:n:p:h",\%opt);
die $usage if ($opt{h});



my $id = $opt{n};
my $pre = $opt{p};
my $sam = $opt{s};


open my $fh,">$pre.$id.informative.sam" or die $!;

my %reads;      # hash to assistant to detect paired reads in sam file
my @rs;			# array contain the pairs of one fragment

open SAM, $sam or die $!;

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
		print_clu ( @rs) if (eof(SAM));
	}
}
sub print_clu{
	my @va = @_;
	my $pt = join "\n",@va;     # $pt is ready to print
	#my $map_q;     # used to check map_q of reads at genome
	
	my $cross;      # used to check if pairs are located at genome and te
	foreach my $it (@va){
		my ($title,$flag,$chr,$mq,$cig,$rnext) = (split /\t/,$it)[0,1,2,4,5,6];
		print "ERROR:$_\n$title,$flag,$chr,$cig,$rnext\n" unless($rnext);    # used to help check the vality of code

		if(($rnext =~ /$id/ and $chr !~ /$id/) or ($rnext ne "=" and $chr =~ /($id)/)){
			$cross ++;
		}
			
		#if ($chr !~ /$id/ and $mq >0){
		#	$map_q ++ ;
		#}
	}		
	
	if($cross){
		print $fh "$pt\n";
	}
}

