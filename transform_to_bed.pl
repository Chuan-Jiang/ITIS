#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

############## parameters ##############
my %opt;
my $usage = "$0
	-p : project_name/ will be prefix of your output file
	-i : <REQUIRED> reads list file support insertion  
	-w : default 100, windows size to cluster reads
	
FILTER parameters:
	-t : defualt [3,1,1], the required number of events support insert: total 3, at TE start 1, ant TE end 1; 
		 If Sites with total larger than or equal 3, it will be saved as raw list \$proj.raw.inser.bed
	-b : original bam file; used to calculate the bg depth
	-r : default 0.2; the ratio of number of support reads with background depth aroung 200bp, assume '-b'
	-D : default <3,200>; the required depth region; assume '-b '
	
	-h : help
	";

die $usage if (@ARGV == 0);
getopts("p:i:t:r:D:b:w:h",\%opt);
die $usage if ($opt{h});

my $proj = $opt{p};

my $ins_file = $opt{i};

my $window = $opt{w}?$opt{w}:100;

my ($total_reads,$te_s,$te_e) = $opt{t}?(split ',',$opt{t}):(split ',','3,1,1');

my $bam = $opt{b}?$opt{b}:0;

my $ratio = $opt{r}?$opt{r}:0.2;

my ($min_d,$max_d) = $opt{D}?(split ',',$opt{D}):(split ',',"3,200");
###############parameters################
#########################################
#########################################



open INS,$ins_file or die $!;
open OUT_F, ">$proj.filtered.bed" or die $!;
open OUT_R, ">$proj.raw.bed" or die $!; 

chomp (my @sites = <INS>);   # put all sorted ins_file in in array

my @list;    # contain all the support reads  list without filtered

for (my $i = 0;;){
	my @clu;
	my $step ;
	push @clu,$sites[$i];        #  $i the the start point of one string of pos
	
	###############  detect a strint of pos at genome #################
	
	my($id_o,$dir_o,$chr_o,$pos_o,$ty_o) = split /\t/,$sites[$i];
	for (my $j =1;;$j++){               # step in to  the next ins pos
		$step = $j;
		my $pre = $i + $j -1;                  
		my $nex = $i + $j;
		last if ($nex >= @sites);	
		my ($id_p,$dir_p,$chr_p,$pos_p,$ty_p) = split /\t/,$sites[$pre];
		my ($id_n,$dir_n,$chr_n,$pos_n,$ty_n) = split /\t/,$sites[$nex];
		
	
		
		my $win_s;  # determine using which window step
		if ($ty_p =~ /C/ and $ty_n =~ /C/){
			$win_s = $window ;
		}elsif($ty_p =~ /G|T/ and $ty_n =~ /G|T/){
			$win_s = 20;
		}else{
			$win_s = $window;
		}
				
		
		if (($dir_p eq $dir_n or $dir_n eq "NA" or $dir_p eq "NA") and  ($chr_p eq $chr_n ) and  (($pos_n - $pos_p ) <= $win_s)){
			push @clu,$sites[$nex];       # @clu have a cluster of support reads each within a window 50bp
		}else{
			last;
		}
	}

	############### use subfunction to determine if there is a TE insertion site ###
	bed(@clu);
	###############   end   #########################

	
	$i = $i + $step;
	last if ($i >= @sites);
}



sub bed{          # use a cluster of support reads to determine if it is a ture TE insertion
	
	my @clu = @_;
	my $sc = scalar @clu;      # $sc    :     the number of support reads
	my %sites;   # contain all start and end site

	my @sit_s;   # start site
	my @sit_e;   # end site
	my @rou_s;   # uncertain start site
	my @rou_e;	 # uncertain end site
	
	my %in_ha;
	my ($id,$dir,$chr,$pos,$ty);
	foreach my $pos (@clu){
		($id,$dir,$chr,$pos,$ty) = split /\t/,$pos;
		$in_ha{$id} ++;
		if ($ty =~ /GS|TS/ ){
			push @sit_s,$pos;
		}elsif ($ty =~ /GE|TE/){
			push @sit_e,$pos;
		}elsif( $ty =~ /CS/){
			push @rou_s,$pos;
		}elsif( $ty =~ /CE/){
			push @rou_e,$pos;
		}
	}
	my $num_fg = keys %in_ha;

	############### collet exact insert pos 
	my $ss;	
	if (@sit_s){
		$ss = "+".mode(@sit_s);
	}elsif(@rou_s) {
		$ss = "-".median(@rou_s);
	}else{
		$ss = 0;
	}
	
	my $ee;
	if (@sit_e){
		$ee = "+".mode(@sit_e);
	}elsif (@rou_e){
		$ee = "-".median(@rou_e);
	}else{
		$ee = 0;
	}	
	##################
	my($clp_s,$clp_e,$crs_s,$crs_e) = (scalar@sit_s,scalar@sit_e,scalar@rou_s,scalar@rou_e);
	my $total = $clp_s+$clp_e+$crs_s+$crs_e;
	
	my $boo = 0;

	if ($total >=  $total_reads and  ($clp_s + $crs_s) >= $te_s and ($clp_e+$crs_e) >= $te_e){     #####i########################  filtering 1#####  filtering 1
		$boo = 1;
	}
	
	my $t = join ",",$clp_s,$clp_e,$crs_s,$crs_e;  #  in the order of 'Reads support Start and End'. Fragment suported Start and End
		
	my($s_p,$e_p);	 #  to print 
	if ($dir eq "R"){
		($s_p,$e_p) = deter_ord($ss,$ee);
		$t = "-$sc=$t";
	}elsif($dir eq "S"){
		($s_p,$e_p) = deter_ord($ee,$ss);
		$t = "+$sc=$t";
	}


	if ($bam){	
		my $pad = int((100 - ($e_p-$s_p))/2);
	
		# $s_r and $e_r are a large region to count the fragment cover the putative insertion site
		# $s_p and $e_p are putative insertion sites
		my ($s_r,$e_r)=$pad>0 ? ($s_p-$pad,$e_p+$pad):($s_p,$e_p);
	
		######## pick the bg depth ######
		#system("samtools view $bam  $chr:$s_r-$e_r -b -h | samtools sort -n - tmp_for_itit.sort 2>/dev/null")== 0 or die " PIPE error when calculate bg reads\n";
	
		open SV,"samtools depth  $bam  -r $chr:$s_r-$e_r  2>/dev/null| " or die $!;
		my $dep;
		my $m_l;
		while(<SV>){
			my ($id,$pos,$d)= (split /\t/,$_)[0,1,2];
			$dep += $d;
			$m_l++;
		}
		if( $m_l and $dep){
			$dep = int($dep/$m_l);
		}else{
			$dep = 0;
		}
		$t .= ";$num_fg/$dep";
		if ( $boo and $dep >= $min_d  and $dep <= $max_d and $num_fg/$dep >= $ratio ){                               ############################Filtering 2#################
			$boo = 1;
		}else{
			$boo = 0;
		}
	}
	
	print OUT_R "$chr\t$s_p\t$e_p\t$t\n" if ($total >= $total_reads);
	print OUT_F "$chr\t$s_p\t$e_p\t$t\n" if ($boo);
}

sub deter_ord{
	my ($aa,$bb) = @_;
	my ($ar,$br);
	if ($aa * $bb == 0){     # if juction of one side is determined
		$ar = ($aa != 0)?abs($aa):abs($bb);
		$br = ($bb != 0)?abs($bb):abs($aa);
		return($ar-1,$br);
	}
  
	if ( abs($aa) <=  abs($bb)){       # if the relation of correct
		($ar,$br) = (abs($aa),abs($bb));
		return($ar-1,$br);
	}
	$ar = $aa>0?$aa:$bb;     # if no correct,
	$br = $bb>0?$bb:$aa;
	($ar,$br) = sort {$a <=> $b} (abs($ar),abs($br));
	return ($ar-1,$br);
}


sub mode {
	my @nums = @_;
	my %ha;
	foreach my $i (@nums){
		$ha{$i} ++;
	}
	my ($mode) = (sort {$ha{$a} <=> $ha{$b}} keys %ha)[-1];
	return $mode;
}

sub median {
	my @num = @_;
	@num = sort {$a <=> $b} @num;
	my $len = scalar @num;
	my $me = $num[int($len/2)];
	return $me;
}



