#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

my %opt;
my $usage =  "USAGE:
	$0 
	-s sam file of short reads aligned to reference  if not specifed , it will read from stdin 
	-n the ID of TE seq 
	-o out_put file prefix
	-h print this help
	" ;


die $usage if ( @ARGV == 0);
getopts("s:n:o:h",\%opt);
die $usage if ($opt{h});



my $id = $opt{n};
my $pre = $opt{o};
my $sam = $opt{s};


open my $fh,">$pre.informative.sam" or die $!;

my %reads;

open SAM, $sam or die $!;

my @rs;	
while (<SAM>){
	chomp;
	if (/^@/){
		print $fh "$_\n";
		next;
	}
	my ($title,$chr,$cig,$rnext) = (split /\t/,$_)[0,2,5,6];
	if (! keys(%reads )  or exists ($reads{$title})){
		push @rs,$_;
		$reads{$title} = "1";
	}else{
		my @va = @rs;
		my $pt = join "\n",@va;
		foreach my $it (@va){
			my ($title,$flag,$chr,$cig,$rnext) = (split /\t/,$it)[0,1,2,5,6];
			print "ERROR:$_\n$title,$flag,$chr,$cig,$rnext\n" unless($rnext);
			if(($rnext =~ /$id/ and $chr !~ /$id/) or ($rnext ne "=" and $chr =~ /($id)/)){
				print $fh "$pt\n";
				last;
			}
		}

		undef(%reads);
		$reads{$title} = 1;
		@rs = ($_);
		#$reads{$title} = "############$_";
	}
}
