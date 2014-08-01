#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;

############## parameters ##############
my %opt;
my $usage = "$0
	-p : project_name/ will be prefix of your output file
	-i : <REQUIRED> reads list file support insertion  
	-w : default 100, windows size to cluster reads
	-n : the name of you te seq	
	-l : the lenght of library insertion
	-b : original bam file; used to calculate the bg depth; bam file should been sorted and indexed
			All reads should been aligned to mergered reference genome file with bwa mem
	-h : help
	";

die $usage if (@ARGV == 0);
getopts("p:i:b:w:n:l:h",\%opt);
die $usage if ($opt{h});

my $proj = $opt{p};

my $ins_file = $opt{i};

my $window = $opt{w}?$opt{w}:100;

my $lib_l = $opt{l}?$opt{l}:500;

my $bam = $opt{b}?$opt{b}:0;

my $te = $opt{n};

###############parameters################
#########################################
#########################################



open INS,$ins_file or die $!;




###################################
####### read through ins.loc.lst###
###################################
###################################

# to get this two variables:
my @lsts;
my @tsds;
# this two variables will be used in the final step

my %rcder;
my @clus;
my %dirs;
while(<INS>){	# iterate the lst file and get tsd information and candidate insertions 
	chomp;
	my ($id,$dir,$chr,$pos,$ty) = split /\t/;

	my $win_s;# determine using which window step

	if(! %rcder ){
		$rcder{chr} = $chr;
		$rcder{pos} = $pos;
		$rcder{ty} = $ty;
		
		$dirs{$dir} ++;
		push @clus,$_;
		next;
	}

	# determine the window step to cluster reads
	if ($ty=~ /C/ and $rcder{ty} =~ /C/){
		$win_s = $window ;
	}elsif($ty =~ /G|T/ and  $rcder{ty} =~ /G|T/){
		$win_s = 20 ;
	}else{
		$win_s = $window;
	}
	# get a cluter of support reads
	if(($chr eq $rcder{chr}) and ($pos - $rcder{pos} <= $win_s)){
		push @clus,$_;
		$rcder{chr} = $chr;
		$rcder{pos} = $pos;
		$rcder{ty} = $ty;
		$dirs{$dir} ++;

		if(eof(INS)){
			my($tsd,$inf) = bed(\@clus,\%dirs);
			push @tsds,$tsd;
			push @lsts, $inf;
		}
	}else{
		my($tsd,$inf) = bed(\@clus,\%dirs);
		push @tsds,$tsd;
		push @lsts, $inf;

		undef(%dirs);
		undef(%rcder) ;
		undef(@clus);
	}
}
close( INS );

######### To the final step to print  each candidate#######

# get the most of tsd length
my $tsd_l = mode(@tsds);
print STDERR "Estimate the length of TSD is 5 bp\n";
$tsd_l--;

# check each  candidate insertion
for my $it ( @lsts){   # iterate to integrate other information

	my($chr,$ss,$ee,$in_ha_ref,$sc,$clp_s,$clp_e,$crs_s,$crs_e,$tags,$dir) = @$it;
	if($dir eq "."){
		print STDERR "Can't determing the direction of insertion at $chr\t$ss\t$sc\n";
		next;
	}
	
	### esitimate the exact insertion site  ####
	my ($s_p,$e_p);
	if($clp_s and $clp_e){
		($s_p,$e_p) = sort {$a<=>$b} ($ss,$ee);
	}elsif($clp_s){
		my $form = $dir.$tsd_l;
		
		($s_p,$e_p) = sort {$a<=>$b} ($ss,$ss-$form);
	}elsif($clp_e){
		my $form = $dir.$tsd_l;
		
		($s_p,$e_p) = sort {$a<=>$b} ($ee,$ee+$form);
	}else{
		($s_p,$e_p) = deter_ord($ss,$ee);
	}

	##  if bam exists, calculate bg depth and esimate homo or hetero
	if ($bam){		
	
		###### esitmate ratio #####
		estimate_homo($chr,$s_p,$e_p,$in_ha_ref) if($clp_s or $clp_e);
		

		######## pick the bg depth ######
		my $pad = int((100- ($e_p-$s_p))/2);
		my ($s_r,$e_r)=$pad>0 ? ($s_p-$pad,$e_p+$pad):($s_p,$e_p);
		
		open SV,"samtools depth  $bam  -r $chr:$s_r-$e_r  2>/dev/null| " or die $!;
		my $dep;
		my $m_l;
		while(<SV>){
			my ($id,$pos,$d)= (split /\t/,$_)[0,1,2];
			$dep += $d;
			$m_l++;
		}
		close SV;
		if( $m_l and $dep){
			$dep = int($dep/$m_l);
		}else{
			$dep = 0;
		}
	}




######################
#######################

sub bed{          # use a cluster of support reads to determine if it is a ture TE insertion

	my @clu  =  @{$_[0]};
	my %dirs =  %{$_[1]};
	
	
	## use ratio to determine the direction of insertion
	my $dir;
	if ($dirs{S} and $dirs{R}){	
		if($dirs{S}/($dirs{R}+$dirs{S}) > 0.8){
			$dir = "+";
		}elsif( $dirs{R}/($dirs{R}+$dirs{S}) > 0.8){
			$dir = "-";
		}else{
			$dir = ".";
		}
	}else{
		($dir) = keys %dirs;
		$dir = ($dir eq "S")?"+":"-";
	}


	###		determine the insertion site
	my $sc = scalar @clu;      # $sc    :     the number of support reads

	my @sit_s;   # exact start site
	my @sit_e;   # exact end site
	my @rou_s;   # uncertain start site
	my @rou_e;	 # uncertain end site
	my @te_s;    # the start site of junction at te
	my @te_e;    # the end site of junction at te

	my %in_ha;   # used to store read id
	my ($id,$d,$chr,$pos,$ty);
	my $map_q;
	foreach my $pos (@clu){
		($id,$d,$chr,$pos,$ty,my $mq) = split /\t/,$pos;
		$in_ha{$id} ++;
		if ($ty =~ /(GS:|TS:)(\d+)/ ){
			push @sit_s,$pos;
			push @te_s, $2;
		}elsif ($ty =~ /(GE:|TE:)(\d+)/){
			push @sit_e,$pos;
			push @te_e, $2;
		}elsif( $ty =~ /CS/){
			push @rou_s,$pos;
		}elsif( $ty =~ /CE/){
			push @rou_e,$pos;
		}
		$map_q += $mq;
	}
	$map_q = int($map_q/@clu);

	my $num_fg = keys %in_ha;    # the total num of read pairs suppor insertion
	             
	### collet exact insert pos 
	
	my $ss;	  # start site at genome
	if (@sit_s){
		$ss = mode(@sit_s);
	}elsif(@rou_s) {
		$ss = median(@rou_s);
	}else{
		$ss = 0;
	}

	my $ee;   #  end site at genome
	if (@sit_e){
		$ee = mode(@sit_e);
	}elsif (@rou_e){
		$ee = median(@rou_e);
	}else{
		$ee = 0;
	}	

	my $t_s;   # start site at TE
	if(@te_s){
		$t_s = mode(@te_s);
	}else{
		$t_s = "NA";
	}
	my $t_e;   # end site at TE
	if(@te_e){
		$t_e = mode (@te_e);
	}else{
		$t_e = "NA";
	}

	##################
	my($clp_s,$clp_e,$crs_s,$crs_e) = (scalar@sit_s,scalar@sit_e,scalar@rou_s,scalar@rou_e);
	my $SR = join ",",$sc,$clp_s,$clp_e,$crs_s,$crs_e;  #  in the order of 'Reads support Start and End'. Fragment suported Start and End
		
	return(abs($ss-$ee),[$chr,$ss,$ee,\%in_ha,$sc,$clp_s,$clp_e,$crs_s,$crs_e,"MQ=$map_q;NM=$te;TS=$t_s;TE=$t_e",$dir]);
}




sub deter_ord{
	my ($aa,$bb) = @_;
	my ($ar,$br);
	if ($aa * $bb == 0){     # if juction of one side is determined
		$ar = ($aa != 0)?abs($aa):abs($bb);
		$br = ($bb != 0)?abs($bb):abs($aa);
		return($ar-1,$br);
	}else{
		($ar,$br) = sort {$a <=> $b} (abs($ar),abs($br));
		return ($ar-1,$br);
	}
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



sub estimate_homo {        # check each read pair  around  the candidate insert sites
	my($chr,$s_r,$e_r,$in_ha_ref) = @_;	
	
	my %in_ha = %$in_ha_ref;
	my $sam_s = $s_r - $lib_l;
	my $sam_e = $e_r + 2;

	#system("samtools view $bam  $chr:$sam_s-$sam_e -b -h | samtools sort -n - $proj/tmp_file_for_homo.sort  2>/dev/null")== 0 or die " PIPE error when calculate bg reads\n";
	#open my $homo "samtools view $proj/tmp_file_for_homo.sort.bam | " or die $!;
	
	###iterate each read pair
	my %reads;
	open my $home "samtools view $bam $chr:$sam_s-$sam_e | " or die $!;
	my @rds;
	while(<$homo>){
		chomp;
		next if (/^@/);
		
		my ($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$_)[0,1,2,3,4,5,6,7,9];

		if (! keys %reads or exists $reads{$id}){
			push @rds , $_;
			if(eof($fh)){
				#####
				check_reads_pair(\@rds);		
			}
		}else{

			######
			check_reads_pair(\@rds);
			######
			
			undef(%reads);
			undef(@rds);
			$reads{$id} = 1;
			push @rds , $_;
		}
	}

}

sub check_reads_pair {
	my @rds = shift @_;
	
	my @range;
	for my $r ( @rds){
		my ( $id,$chr,$pos,$nchr,$npos,$tlen) = (split /\t/,$r)[0,2,3,6,7,8];
		if (exists $in_ha{$id}){
			return "I";
		}else{
			if ($chr =~ /$te/ or $chr =~ /$te/){
				return "I";
			}elsif($chr !~ /$te/ and $nchr eq "="){
				if($tlen > 0){
					@range = ($pos,$pos+$tlen);
				}else($tlen < 0){
					@range = ($npos,$npos-$tlen) ;
				}
			}else{

		}
		return \@range;
	}
}









