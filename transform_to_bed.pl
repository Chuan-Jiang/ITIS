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
	-d : swithh on debug model
	";

die $usage if (@ARGV == 0);
getopts("p:i:b:w:n:l:hd",\%opt);
die $usage if ($opt{h});

my $proj = $opt{p};

my $ins_file = $opt{i};

my $window = $opt{w}?$opt{w}:100;

my $lib_l = $opt{l}?$opt{l}:500;

my $bam = $opt{b}?$opt{b}:0;

my $te = $opt{n};

my $db = $opt{d};
###############parameters################
#########################################
#########################################



open INS,$ins_file or die $!;
open OUT_R, ">$proj.raw.bed" or die $!;



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
		rcd_it($_,$chr,$pos,$ty,$dir);	
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
		rcd_it($_,$chr,$pos,$ty,$dir);

		if(eof(INS)){
			my($tsd,$inf) = collect_infor(\@clus,\%dirs);
			print "CLU:@clus\n" if $db;
			push @tsds,$tsd if $tsd>0;
			push @lsts, $inf;
		}
	}else{
		my($tsd,$inf) = collect_infor(\@clus,\%dirs);
		print "CLU:@clus\n" if $db;
		push @tsds,$tsd if $tsd >0;
		push @lsts, $inf;

		undef(%dirs);
		undef(%rcder) ;
		undef(@clus);
		rcd_it($_,$chr,$pos,$ty,$dir);
		if (eof(INS)){
			print "CLU:@clus\n" if $db;
			my($tsd,$inf) = collect_infor(\@clus,\%dirs);
			push @tsds,$tsd  if $tsd >0;
			push @lsts, $inf;
		}
	}
}
close( INS );

sub rcd_it{
	my($it,$chr,$pos,$ty,$d) = @_;
	push @clus, $it;
	$dirs{$d} ++;
	$rcder{chr} = $chr;
	$rcder{pos} = $pos;
	$rcder{ty} = $ty;
}

######### the final step to print  each candidate#######
#  1, used esimated tsd length to determine the exact insertion site
#
#  2,calculate the bg depth and estimate Genome Type if BAM file is used  #
##########################


# get the most of tsd length
my $tsd_l = mode(@tsds);
$tsd_l++;
print STDERR "Estimate the length of TSD is $tsd_l bp\n";
$tsd_l--;


# check each  candidate insertion
for my $it ( @lsts){   # iterate to integrate other information

	my($chr,$ss,$ee,$in_ha_ref,$sc,$clp_s,$clp_e,$crs_s,$crs_e,$tags,$dir) = @$it;
	if($dir eq "."){
		print STDERR "Can't determing the direction of insertion at $chr\t$ss\t$sc\n";
		next  ;     # if a cluster of reads have different diretion, discarded
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
		my ($sup,$nsup) = ("NA","NA");
		($sup,$nsup) = estimate_homo($chr,$s_p,$e_p,$in_ha_ref,$dir) if($clp_s or $clp_e);
		
		$tags .= ";GT=$sup,$nsup";
		######## pick the bg depth ######
		
		my $pad = int((100- ($e_p-$s_p))/2);
		my ($s_r,$e_r)=$pad>0 ? ($s_p-$pad,$e_p+$pad):($s_p,$e_p);
		
		$s_r = 0 if $s_r < 0;

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
		$tags .= ";DP=$dep";
	}
	print "$chr\t$s_p\t$e_p\t$tags\t.\t$dir\n" if $db;
	print OUT_R  "$chr\t$s_p\t$e_p\t$tags\t.\t$dir\n";
}

######################
#######################





# the aim of this subroutine is to collect the information of candidate sites
sub collect_infor{          # use a cluster of support reads to determine if it is a ture TE insertion

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
		$in_ha{$id} .= $ty;
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
	my $SR = join ",", $num_fg,$sc,$clp_s,$clp_e,$crs_s,$crs_e;  #  in the order of 'Reads support Start and End'. Fragment suported Start and End
	
	my $tsd_l = (@sit_s and @sit_e)?abs($ss-$ee)	: 0;
	return($tsd_l,[$chr,$ss,$ee,\%in_ha,$sc,$clp_s,$clp_e,$crs_s,$crs_e,"SR=$SR;MQ=$map_q;NM=$te;TS=$t_s;TE=$t_e",$dir]);
}

# SR : counts of total and every type of supporting reads
# MQ : the average mapping quality
# NM : the name of te 
# TS : the start site of te 
# TE : the end site of te



sub deter_ord{
	my ($aa,$bb) = @_;
	my ($ar,$br);
	if ($aa * $bb == 0){     # if juction of one side is determined
		$ar = ($aa != 0)?abs($aa):abs($bb);
		$br = ($bb != 0)?abs($bb):abs($aa);
		return($ar-1,$br);
	}else{
		($ar,$br) = sort {$a <=> $b} (abs($aa),abs($bb));
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
	my($chr,$s_r,$e_r,$in_ha_ref,$dir) = @_;	
	print "\nEstimate:$chr\t$s_r\t$e_r\n" if $db;
	my %in_ha = %$in_ha_ref;
	my $sam_s = $s_r - $lib_l;
	my $sam_e = $e_r + $lib_l;

	$sam_s = 0 if $sam_s < 0;	
	
	if($db){
		while(my($k,$v) = each %in_ha){
			print "SUPP: $k\t$v\n";
		}
	}	
	###iterate each read pair
	my %reads;
	print "SAMVIEW:samtools view -X $bam $chr:$sam_s-$sam_e\n" if $db;
	open my $sam, "samtools view -X $bam $chr:$sam_s-$sam_e | " or die $!;
	while(<$sam>){
		chomp;
		next if (/^@/);
		my $r = $_;
		my ( $id,$tag,$chr,$pos,$mq,$cig,$nchr,$npos,$tlen) = (split /\t/,$r)[0,1,2,3,4,5,6,7,8];
		print "RD: $r\n" if $db;
		#next if ($mq == 0);
		
		next if ($tag =~ /u/);
		$reads{$id} = 9;

		#  1 :  support insertion
		#  2 :  support excision
		#  3 :  flank reads 
		#  4 :  discarded support reads

		if (exists $in_ha{$id} and $reads{$id} >= 1){
			 $reads{$id} = 1;
		}else{
			if($nchr eq "=" and $tlen != 0 and abs($tlen) < 3*$lib_l ){	
				
				
				my @range;	#  determine the range of pair
				
				if($tlen > 0){
						
					@range = ($pos,$pos+$tlen-1);
				}elsif($tlen < 0){
					@range = ($npos,$npos-$tlen-1) ;
				}

				# check if the range overlap with TE insertion site
				if ($range[0] < $s_r and $range[1] > $e_r and $reads{$id} >= 2 ){
					$reads{$id} = 2;
				}else{
					$reads{$id} = 4 if ($reads{$id} >= 4);
				}
			}elsif($nchr eq $te){
				$reads{$id} = 3 if ($reads{$id} >= 3);
			}else{
				$reads{$id} = 4 if ($reads{$id} >= 4);
			}
		}
	}
	my $sup = 0;
	my $nsup = 0 ;
	my (@sup_r,@nsup_r) if $db;
	while(my ($k,$v) = each %reads){
		if($v == 1){
			push @sup_r, $k if $db;
			$sup++;
		}elsif($v == 2){
			push @nsup_r, $k if $db;
			$nsup ++;
		}
	}
	if($db){
		print "SUP: @sup_r\n";
		print "NSUP: @nsup_r\n";
	}
	return ($sup,$nsup);

}









