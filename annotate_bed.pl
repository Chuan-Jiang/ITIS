#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Getopt::Std;

my %opt;

my $usage =   "USAGE:
	$0 :
	-b bed file of insertion locs {REQUIRED}
	-a annotation file in gff3 format { If not provied , only generate snpbatch filt for IGV}
	-g merged genome sequence with TE          {REQUIRED}
	-n the name of your retrotransposon name   {REQUIRED}
	-p the prefix you project                  {REQUIRED}
	-h print this help infor                    
	-d the directory contain your results      Default: Current directory 
	" ;
die $usage if (@ARGV == 0);
getopts("b:a:g:i:n:p:d:h",\%opt);
die $usage if ($opt{h});
my($bed_file,$gff_file,$genome_file,$te,$proj) = ($opt{b},$opt{a},$opt{g},$opt{n},$opt{p});
my $folder = $opt{d}?$opt{d}:'.';
my $wd = $ENV{'PWD'};

################ generate igv batch   #############
	open FH,"> $folder/$proj.$te.igv.bat" or die $!;
	open BED,"$bed_file" or die $!;
	print FH "snapshotDirectory $wd/$proj.$te.snap.dir\n";
	while (<BED>){
		chomp;
		 my ($chr,$s,$e) = split /\t/,$_;
		 if (($e - $s) > 300){
			 print STDERR "A large region: $chr,$s,$e\n";
			 print FH "goto $chr:$s-$e\n";
			 print FH "snapshot ${chr}_${s}_${e}.png\n";
		 }else{
			 my $slop = 150 - ($e - $s);
			 my $s_p = $s - int($slop/2);
			 my $e_p = $e + int($slop/2);
			 print FH "goto $chr:$s_p-$e_p\n";
			 print FH "snapshot ${chr}_${s}_${e}_slop$slop.png\n";
		 }
	 }

if($gff_file){    #  if provided gff file, then process the following code util the end

################# generate usefull gff file containg intergenic region  ##########

open OUT, ">$folder/$proj.$te.annotation.tsv" or die $!;
open TEM,">tem.$proj" or die $!;
open GFF,"awk 'BEGIN{IGNORECASE=1} {if(\$3 ~ /gene/){print \$0}}' $gff_file |" or die $!;
chomp (my @gff = <GFF>);
for (my $i = 1;$i<@gff;$i++){
	print TEM "$gff[$i-1]\n";
	my $l = $gff[$i-1];
	my $n  = $gff[$i];
	my @las = split /\t/,$l;
	my @nos = split /\t/,$n;
	if (($nos[3] - $las[4]) > 0 and $nos[0] eq $las[0]){
		my $s = $las[4] +1;
		my $e = $nos[3] -1;
	
		my %la;
		foreach  (split /;/,$las[8]){
			my ($k,$v) = split /=/,$_;
			$la{$k} = $v;
		}

		my %ne;
		foreach (split /;/,$nos[8]){
			my ($k,$v) = split /=/,$_;
			$ne{$k} = $v;
		}

		my $la_ge = $la{ID};
		my $la_no = $la{Note};

		my $ne_ge = $ne{ID};
		my $ne_no = $ne{Note};
		print TEM "$las[0]\t$las[1]\tIntergenic\t$s\t$e\t.\t.\t.\tID=${la_ge}_${ne_ge};Note=Intergenic_${la_no}_INSERT_$ne_no\n";
	}
}
##################### put genome in hash##################################

my %genome = Seq::seq_hash($genome_file);
##########################################################################



################### intersect using bedtools  ###########################

open BEDTOOL, "bedtools intersect -a $bed_file -b tem.$proj -wa -wb |" or die $!;

my $la=0;
my $la_t;
while (<BEDTOOL>){
	chomp;
	#print "$_\n";
	my($chr,$s,$e,$name,$d,$t,$an) = (split /\t/,$_)[0,1,2,3,5,8,14];
	
	my ($up,$down);

	my $te_seq = ($d =~ /\+/)?$genome{$te}:Seq::rev_com($genome{$te});
	
	$up = substr ($genome{$chr},$e-1005,1005);
	$down = substr ($genome{$chr},$s,1005);
	
	my $seq = $up.lc($te_seq).$down;
	
	print OUT "$chr\t$s\t$e\t$name\t$d\t$t\t$an\t$seq\n";

	$la = "$chr:$s:$e";
	$la_t = $t;

}

#unlink "tem.$proj";

}
