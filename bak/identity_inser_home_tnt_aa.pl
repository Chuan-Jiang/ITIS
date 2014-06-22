#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Bio::Seq;
my ($sam_file,$genome_file,$ins_size,$te_idx,$te,$ltr) = @ARGV;
print "$0 <sam> <genome> <ins_size> <te_idx> <te_name> <ltr_len>\n" unless @ARGV;


############# put genome seq in hash  #########
my %genomes = Seq::seq_hash($genome_file);
my $tnt_len = length ($genomes{$te});  # transposon length


my $tar;
my $rds;
my $num;
my $num_eff;
my $num_f;
my $rex_te = '^(\d+[SH]\d+M|\d+M|\d+M\d+[SH])$';

###########################################################
################ the mainbody of code######################
###########################################################
###########################################################


my @aligns = scan_sam ($sam_file);
foreach my $read (@aligns){
	my @hits = split /############/,$read;
	my ($co,$r_h) = find (@hits );
	next if ($co == 0);
	my %cors = %$co;

	my $tar_cig = $cors{tar}{cig};
	my $te_cig = $cors{$te}{cig}; 
	
	unless ($tar_cig and $te_cig){
		print "can not determine the pairs\n";
		die;
	}

	if ( $tar_cig =~ /M/ and $te_cig =~ /M/){
		cross(%cors);
	}
	if($te_cig =~ /S/){
		te_start(%cors);
	}
	if($te_cig =~ /E/){
		te_end(%cors);
	}
	if($tar_cig =~ /S/){
		ge_start(%cors);
	}
	if($tar_cig =~ /E/){
		ge_end(%cors);
	}
	if ($tar_cig =~ /Z/){
		print "tar_cig_f:$read\n";
	}
	if ($te_cig =~ /Z/){
		print "te_cig_f:$read\n";
	}
}

#################### the mainbody of code ######################


sub scan_sam{
	my $file = shift @_;
	open my $fh , $file or die $!;
	my %reads;
	my @re;
	while (<$fh>){
		chomp;
		next if (/^@/);
		my $id = (split /\t/,$_)[0];
		if(! keys %reads){
			$reads{$id} = "$_";
		}elsif( exists $reads{$id}){
			$reads{$id} .= "############$_";
		}else{
			my ($g) = values %reads;
			push @re, $g;
			undef(%reads);
			$reads{$id} = "$_";
		}
	}
	return @re;
}	


sub find{
	$num ++;
	my @reads = @_;
	$rds = join "\n",@_; ## in order to print informative reads in one file 
	
	my %cors;
		
	################# put reads in hash #########

	my %reads_hash;
	
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9];	
		(my $r = $flag ) =~ s/\w+(\d)s?/$1/;
		my $p = ($chr =~ /$te/)?$te:"chr";
		$reads_hash{$r}{$p}{pos} = $pos;
		$reads_hash{$r}{$p}{chr} = $chr;
		$reads_hash{$r}{$p}{cig} = $cig;
	}

	##############determine the reads at tnt1#############
	my %trans;  # select the most promising tnt reads
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9];
		#print "$id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq\n";
		(my $r = $flag ) =~ s/\w+(\d)s?/$1/;
		my $pair = ($r == 1)?2:1;
		my $rc = (($flag =~ /r/)? -1:1); 
		
		if (($chr =~ /$te/) and exists ($reads_hash{$pair}{chr} )){
			
			my $num_of_m;
			while ($cig =~ /(\d+)M/g){
				$num_of_m += $1;
			}
			if (exists $trans{nm} and $trans{nm} > $num_of_m){
				next;
			}
			
			$trans{chr} = $chr;
			$trans{nm} = $num_of_m;
			$trans{nchr} = $reads_hash{$pair}{chr}{chr};
			$trans{pos}  = $pos;
			$trans{npos} = $reads_hash{$pair}{chr}{pos};
			$trans{ori} = (($cig =~ /H/)?$r*$rc:0);
			
			if ($cig =~ /H/){
				$seq = comp($trans{ori},\@reads);
			}
			$cors{$te}{seq} = $seq;
			$cors{$te}{id} = $id;
			$cors{$te}{direc} = $rc;
		}	
	}
	
	######### save  pair alignment ###########
	############################################
	############################################
	
	

	################  if tnt have long ltr execute following code#############
	my %alns;   ## key will be position;
	open BT, ("bowtie2 -x $te_idx -c '$cors{$te}{seq}' --local -a --quiet |samtools view -S -X - 2>/tmp/sam.tmp|") or die $!;
	my $aln_sc = 10;
	##print "ididid :: $cors{$te}{id}\n";
	while (<BT>){
		next if (/^@/);
		chomp;	
		my ($pos,$cig,$as) = (split /\t/,$_)[3,5,11];
		my $s = $1 if ($as =~ /AS:i:(\d+)/);
		unless ($s){
			return (0,0);
		}	
		#	print "aln:$_\n";
		if ($s > ($aln_sc - 10)){
			$alns{$pos} = $cig;
		}
		$aln_sc = $s;
	}

		
	my $te_pos;
	if ($cors{$te}{direc} == 1){
		($te_pos) = (sort {$a<=>$b} keys %alns)[-1];
	}else{
		($te_pos) = (sort {$a<=>$b} keys %alns)[0];
	}
	my $te_cig = $alns{$te_pos};
	my $te_cs = cigar($te_cig);
	$cors{$te}{pos} = $te_pos;
	$cors{$te}{cig} = $te_cs;


	############## collect information of reads at genome###############
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9]; 
		my $rc = (($flag =~ /r/)? -1:1);    # -1 antisense   and 1 means sense strand
		(my $r = $flag ) =~ s/\w+(\d)/$1/;
		my $cs = cigar($cig);
		#### select informative reads######
		if ($chr eq $trans{nchr} and ($pos == $trans{npos})){
			my $ori = (($cig =~ /H/)?$r*$rc:0);
			if ($cig =~ /H/){
				$seq = comp($ori,\@reads);
			}
			$cors{tar}{cig} = $cs;       # simplified cigar: M S E Z
			$cors{tar}{direc} = $rc;      # reads direction
			$cors{tar}{id} = $id;	       # read id 
			$cors{tar}{pos} = $pos;       ### when reads if reverse complement the start loc of reads should plus the length of reads
			$cors{tar}{seq} = $seq ;
			$cors{tar}{chr} = $chr;
		}
	}
	return (\%cors,\%reads_hash);
}

sub comp{
	my ($ori,$ref) = @_;
	my @reads = @$ref;
	#### if are hard clipped, exacute following code to get the full sequence  at tnt#######
	my $va;
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9];	
		(my $r = $flag ) =~ s/\w+(\d)/$1/;  #  which reads   1 or 2 ?????
		my $rc = (($flag =~ /r/)? -1:1);    # -1 antisense   and 1 means sense strand

		my $test_n = (($cig !~ /H/)?$r*$rc:0);
		if (abs($ori)== abs ($test_n) and ($ori)*($test_n) > 0){
			$va = $seq;
		}elsif(abs($ori)== abs ($test_n)  and ($ori)*($test_n) < 0){
			$va = Seq::rev_com($seq);
		}
	}
	return $va;
}


sub cigar {
	my $cig = shift @_;

	my $len;		        #  reads length
	my $cs;                      #  simplified cigars 
	my ($m_l);					# match length
		
	while ($cig =~ /(\d+)([MSIH])/g){
		my $n = $1;
		my $c = $2;
		$len += $n;
		$m_l += $1 if ($c eq "M");
	}

	if ($m_l/$len > 0.95 ){
		$cs = "M";
	}else{
		if ( $cig =~ /^(\d+)[SH]/ and $1 >= 3){
			$cs .= "S:$1";
		}
		if ( $cig =~ /(\d+)[SH]$/ and $1 >= 3){
			$cs .= "E:$1";
		}else{
			$cs = "Z";
		}
	}
	return $cs;
}
sub cross {
	my %cors = @_;
	if ( $cors{$te}{direc} == 1 and $cors{$te}{pos} > ($tnt_len - $ins_size)){
		my ($ins_direc,$jun);
		if ($cors{tar}{direc} == 1){
			$ins_direc = "R";
			$jun = $cors{tar}{pos} + length ($cors{tar}{seq});
		}else{
			$ins_direc = "S";
			$jun = $cors{tar}{pos};
		}
		print "$cors{$te}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tCE\n";
		print STDERR "$rds\n";
		$num_eff ++;
	}elsif ( $cors{$te}{direc} == -1  and $cors{$te}{pos} < ($ins_size - length($cors{$te}{seq}))){
		my ($ins_direc,$jun);
		if ( $cors{tar}{direc} == 1 ){
			$ins_direc = "S";
			$jun = $cors{tar}{pos} + length($cors{tar}{seq});
		}else{
			$ins_direc = "R";
			$jun = $cors{tar}{pos};
		}
		print "$cors{$te}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tCS\n";
		print STDERR "$rds\n";
		$num_eff++;
	}else{
		print "cross fail:$cors{$te}{id}\tte_direc:$cors{$te}{direc}\tte_pos:$cors{$te}{pos}\ttar_direc:$cors{tar}{direc}\n"
	}
}

sub te_start{
	my %cors = @_;
	if($cors{$te}{cig} =~ /S:(\d+)/ and $cors{$te}{direc} == -1 ){
		my $l = $1;
		if ($cors{$te}{pos} <=2 ){
			my $que = substr($cors{$te}{seq},0,$l);
			( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
			$que_r = reverse $que_r;
			my $ins_direc;
			my $sub;
			my $chr_t = $cors{tar}{chr};
			my $pos_t = $cors{tar}{pos};
			my $len_t = length $cors{tar}{seq};
			if($cors{tar}{direc} == 1){
				$ins_direc = "S";
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_size);
				
				my ($jun,$diff) = mat($que,$sub);
				
				if($jun != -1 and $diff/$l < 0.05){
					$jun = $jun + $pos_t;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTS\n";
					print STDERR "$rds\n";
				}
			}else{
				$ins_direc = "R";
				my $sub_start = $pos_t-$ins_size +$len_t;
				$sub = substr($genomes{$chr_t},$sub_start-1,$ins_size);
				
				my($jun,$diff) = mat($que_r,$sub);

				if($jun != -1 and $diff/$l < 0.05){
					$jun = $jun + $sub_start;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTS\n";
					print STDERR "$rds\n";
				}	
			}
		}else{
			print "te_start_1:$cors{$te}{id}\t$cors{$te}{pos}\n";
		}
	}else{
		print "te_start_faili:$cors{$te}{id}\t$cors{$te}{cig}\t$cors{$te}{direc}\n";
	}
}

sub te_end{
	my %cors = @_;
	if ($cors{$te}{cig} =~ /E:(\d+)/ and $cors{$te}{direc} == 1){
		my $l = $1;
		my $end_pos = $cors{$te}{pos} + length($cors{$te}{seq}) - $l;
		if ($end_pos >= $tnt_len-10 ){
			my $que = substr($cors{$te}{seq},-$l);
			( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
			$que_r = reverse $que_r;
			my $ins_direc;
			my $sub;
			my $chr_t = $cors{tar}{chr};
			my $pos_t = $cors{tar}{pos};
			my $len_t = length $cors{tar}{seq};
			if($cors{tar}{direc} == 1){
				$ins_direc = "R";
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_size);
				
				my ($jun,$diff) = mat($que_r,$sub);
				
				if( $jun != -1 and $diff/$l <0.05 ){
					$jun = $jun + $pos_t;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTE\n";
					print STDERR "$rds\n";
				}else{
					$num_f ++;
				}
			}else{
				$ins_direc = "S";
				my $sub_start = $pos_t+$len_t-$ins_size;
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_size);
				
				my ($jun,$diff) = mat($que,$sub);

				if($jun != -1 and $diff/$l < 0.05 ){
					$jun = $jun + $pos_t;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTE\n";
					print STDERR "$rds\n";
				}else{
					$num_f ++;
				}
			}
		}else{
			print "te_end_1:$cors{$te}{id}\t$end_pos\n";
		}
	}else{
		print "te_end_fail:$cors{$te}{id}\t$cors{$te}{cig}\t$cors{$te}{direc} \n";
	}
}



sub ge_start{
	my %cors = @_;
	if($cors{tar}{cig} =~ /S:(\d+)/){
		my $l = $1;
		my $que = substr($cors{tar}{seq},0,$l);
		( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
		$que_r = reverse $que_r; 
		my $sub_h = "NNNNN".substr($genomes{$te},0,$l+5);
		my $sub_t = substr($genomes{$te},-($l+5))."NNNNN";
		
		my $ma = match($que,$que_r,$sub_h,$sub_t);
		if ($ma != 0){
			if (abs($ma) == 3){
				my $ins_direc = "S";
				my $jun  = $cors{tar}{pos};
				print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGS\n";
				print STDERR "$rds\n";
				$num_eff++;
			}elsif(abs ( $ma )== 2 ) {
				my $ins_direc = "R";
				my $jun = $cors{tar}{pos};
				print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGS\n";
				print STDERR "$rds\n";
				$num_eff++;
			}
		}else{
			print "te_ma:$cors{tar}{id}\t$ma:$que\t$que_r\n";
		}
	}else{
		print "ge_start:$cors{tar}{id}\t$cors{tar}{cig}\n";
	}
}

sub ge_end{
	my %cors = @_;
	if($cors{tar}{cig} =~ /E:(\d+)/){
		my $l = $1;
				
		my $que = substr($cors{tar}{seq},-$l);
		( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
		$que_r = reverse $que_r; 
		my $sub_h = "AAAAA".substr($genomes{$te},0,$l+5);
		my $sub_t = substr($genomes{$te},-($l+5))."AAAAA";
		
		my $ma = match($que,$que_r,$sub_h,$sub_t);
		if ($ma != 0){
			if (abs($ma) == 3){
				my $ins_direc = "R";
				my $jun = $cors{tar}{pos} + length($cors{tar}{seq})-$l;
				print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGE\n";
				print STDERR "$rds\n";
				$num_eff++;
			}elsif ( abs($ma) == 2 ){
				my $ins_direc = "S";
				my $jun = $cors{tar}{pos} + length($cors{tar}{seq})-$l;
				print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGE\n";
				print STDERR "$rds\n";
				$num_eff++;
			}
		}else{
			print "ge_end:$cors{tar}{id}\t$que\t$que_r\n";
		}
	}else{
		print "ge_end:$cors{tar}{id}\t$cors{tar}{cig}\n";
	}
}


sub mat { 
	my ($que,$sub) = @_;
	my $q_l = length $que;
	my $s_l = length $sub;
	my $record = 100;
	my $loc = -1;
	for my $i (0..($s_l-$q_l)){
		my $tgt = substr($sub,$i,$q_l);
		my $diffcount = () = ($que ^ $tgt) =~ /[^\x00]/g;
		if ($diffcount < $record){
			$loc = $i;
			$record = $diffcount;
		}

	}
	return ($record,$loc);
}


sub match {
	my %relas;
	for my $j (0..1){
		for my $k (2..3){
			for my $i (0..10){
				my $que = $_[$j];
				my $q_l = length $que;
				my $sub = substr($_[$k],$i,$q_l);
				my $diffcount = () = ($que ^ $sub) =~ /[^\x00]/g;
				my $ratio = $diffcount/$q_l;
				#print "$diffcount\t$len_q\n";
				my $direc = ($j==0?1:-1);
				$relas{$ratio} = $direc*$k;
			}
		}
	}
	foreach (sort {$a<=>$b} keys %relas){
		if ( $_ <= 0.05){
			my $v = $relas{$_};
			return "$v";   #### return value SH  ST RH  RT  
		}else{
			return "0";
		}
		last;
	}
}

