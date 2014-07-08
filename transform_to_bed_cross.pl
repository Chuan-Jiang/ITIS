#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

############## parameters ##############
my %opt;
my $usage = "$0
	-b : original bam file; used to calculate the bg depth
	-i : reads file support insertion
	-t : defualt 3, the required number of reads support insert
	-w : default 100, windows size to cluster reads
	-h : help

	";

die $usage if (@ARGV == 0);
getopts("b:i:t:w:h",\%opt);
die $usage if ($opt{h});

my $tre = $opt{t}?$opt{t}:3;
my ($ins_file,$bam) = ($opt{i},$opt{b});
my $window = $opt{w}?$opt{w}:100;
###############parameters################
#########################################
#########################################



open INS,$ins_file or die $!;
chomp (my @sites = <INS>);   # put all sorted ins_file in in array

my @list;    # contain all the support reads  list without filtered

for (my $i = 0;;){
	my @clu;
	my $step ;
	push @clu,$sites[$i];        #  $i the the start point of one string of pos

	###############  detect a strint of pos at genome #################
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
		
		
		if (($dir_p eq $dir_n ) and  ($chr_p eq $chr_n ) and  (($pos_n - $pos_p ) <= $win_s)){
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
	my $t = join ",",scalar@sit_s,scalar@sit_e,scalar@rou_s,scalar@rou_e;  #  in the order of 'Reads support Start and End'. Fragment suported Start and End
		
	my($s_p,$e_p);	 #  to print 
	if ($dir eq "R"){
		($s_p,$e_p) = deter_ord($ss,$ee);
		$t = "-$sc=$t";
	}elsif($dir eq "S"){
		($s_p,$e_p) = deter_ord($ee,$ss);
		$t = "+$sc=$t";
	}

	
	my $pad = int((100 - ($e_p-$s_p))/2);
	
	# $s_r and $e_r are a large region to count the fragment cover the putative insertion site
	# $s_p and $e_p are putative insertion sites
	my ($s_r,$e_r)=$pad>0 ? ($s_p-$pad,$e_p+$pad):($s_p,$e_p);

	######## pick the bg depth ######
	#system("samtools view $bam  $chr:$s_r-$e_r -b -h | samtools sort -n - tmp_for_itit.sort 2>/dev/null")== 0 or die " PIPE error when calculate bg reads\n";
	
	open SV,"samtools mpileup $bam  -r $chr:$s_r-$e_r  2>/dev/null| " or die $!;
	my %pos_h;  # save each reads position
	my $dep;
	while(<SV>){
		my ($id,$pos,$cig,$seq)= (split /\t/,$_)[0,3,5,9];
		my($s,$e);
		if ($cig =~ /^(\d+)(S|H)/){
			$s = $pos-$1;
			if($2 eq "S"){
				$e = $s;
			}elsif($2 eq "H"){
				$e = $pos;
			}
		}else{
			$s = $pos;
			$e = $pos;
		}
		my $seq_len = length($seq);
		$e = $e + $seq_len -1;
		if($cig =~ /(\d+)H$/){
			$e = $e+$1;
		}
		if(exists $pos_h{$id}){
			push @{$pos_h{$id}},($s,$e);
		}else{
			$pos_h{$id} = [$s,$e];
		}
	}

	my $num_bg ; # the number of reads at bg
	while (my($k,$v) = each %pos_h){
		my ($h,$t) = (sort {$a<=>$b} @$v)[0,-1];
		if ($h < $e_p and $t > $s_p){
			$num_bg++;
		}
	}

	print "$chr\t$s_p\t$e_p\t$t;$num_fg/$num_bg\n" if ( $sc > $tre and $num_fg/$num_bg );

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



