#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Bio::Seq;
my ($sam_file,$genome_file,$ins_len,$te,$ltr) = @ARGV;
print "$0 <sam> <genome> <ins_size> <te_name> <ltr_len>\n" unless @ARGV;


############# put genome seq in hash  #########
my %genomes = Seq::seq_hash($genome_file);
my $tnt_len = length ($genomes{$te});  # transposon length


my $tar;
my $rds;
my $num;
my $num_eff;
my $num_f;


################ the mainbody of code######################
my @aligns = scan_sam ($sam_file);
foreach my $read (@aligns){
	my @hits = split /############/,$read;
	my %cors = find (@hits );
	
	######### check if alignment at tnt can 
	my @tar_cig = @{$cors{tar}{cig}};
	my @te_cig = @{$cors{$te}{cig}}; 


	print "@tar_cig\n";
	print "@te_cig\n";
	if ( @tar_cig ~~ /M/ and @te_cig ~~ /M/){
		cross(%cors);
	}
	if(@te_cig ~~ /S/){
		te_start(%cors);
	}
	if(@te_cig ~~ /E/){
		te_end(%cors);
	}
	if(@tar_cig ~~ /S/){
		ge_start(%cors);
	}
	if(@tar_cig ~~ /E/){
		ge_end(%cors);
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
	my $cors_ref = \%cors;
	
	my %trans;  # select the most promising tnt reads
	
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9];
		#print "$id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq\n";
		if ($chr =~ /$te/ and $nchr ne "="){
			
			my $num_of_m;
			while ($cig =~ /(\d+)M/g){
				$num_of_m += $1;
			}
			if (exists $trans{nm} and $trans{nm} > $num_of_m){
				next;
			}
			
			$trans{chr} = $chr;
			$trans{nm} = $num_of_m;
			$trans{nchr} = $nchr;
			$trans{pos}  = $pos;
			$trans{npos} = $npos;
			
			(my $r = $flag ) =~ s/\w+(\d)s?/$1/;
			my $rc = (($flag =~ /r/)? -1:1); 
			$trans{cig} = (($cig =~ /H/)?$r*$rc:0);
		}	
	}
	
	######### save  pair alignment ###########
	############################################
	############################################
	
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9];	
		(my $r = $flag ) =~ s/\w+(\d)/$1/;  #  which reads   1 or 2 ?????
		my $rc = (($flag =~ /r/)? -1:1);    # -1 antisense   and 1 means sense strand

		#### if tnt are hard clipped, exacute following code to get the full sequence
		my $test_n = (($cig !~ /H/)?$r*$rc:0);
		if (abs($trans{cig})== abs ($test_n) and ($trans{cig})*($test_n) > 0){
			$cors{$te}{seq} = $seq;
		}elsif(abs($trans{cig})== abs ($test_n)  and ($trans{cig})*($test_n) < 0){
			$cors{$te}{seq} = Seq::rev_com($seq);
		}
			

        ####### select informative reads######
		my $wei = ($chr =~ /$te/)?$te:"tar";
		if ($chr =~ /$trans{nchr}|$trans{chr}/ and ($pos == $trans{npos} or $pos == $trans{pos})){
				
			my $len;		        #  reads length
			my @cs;                      #  simplified cigars 
			my ($m_l);					# match length
			
			while ($cig =~ /(\d+)([MSIH])/g){
				my $n = $1;
				my $c = $2;
				$len += $n;
				$m_l += $1 if ($c eq "M");
			}

			if ($m_l/$len > 0.95 ){
				push @cs, "M";
			}
			if ( $cig =~ /^(\d+)[SH]/ and $1 >= 3){
				push @cs, "S:$1";
			}
			if ( $cig =~ /(\d+)[SH]$/ and $1 >= 3){
				push @cs, "E:$1";
			}
			#print "$id\t$s\n";	
			$cors{$wei}{cig} = \@cs;      # simplified cigar: M S E Z
			$cors{$wei}{direc} = $rc;   # reads direction
			$cors{$wei}{chr} = $chr;    # which chromosome
			$cors{$wei}{id} = $id;	  # read id 
			$cors{$wei}{pos} = $pos;  ### when reads if reverse complement the start loc of reads should plus the length of reads
			$cors{$wei}{seq} = $seq;
			$cors{$wei}{mq} = $mq;
		}
	}
	return %cors;
}


sub cross {
	my %cors = @_;
	if ( $cors{$te}{direc} == 1 and ( $cors{$te}{pos} > ($tnt_len - $ltr) or $cors{$te}{pos} < ($ltr - 100))){
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
	}elsif ( $cors{$te}{direc} == -1 and ($cors{$te}{pos} <  ($ltr-100)  or $cors{$te}{pos} > ($tnt_len - $ltr) )){
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
		$num_f ++;
	}
}

sub te_start{
	my %cors = @_;
	if($cors{$te}{cig} =~ /S:(\d+)/ ){
		my $l = $1;
		if ($cors{$te}{pos} <=2 or ($cors{$te}{pos} <= ($tnt_len - $ltr + 2) and $cors{$te}{pos} > ($tnt_len - $ltr))){
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
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_len);
				my $jun = -1;
				my $m_r = 5 ;
				for my $i (0..(length ($sub)- length ($que))){
					my $s = substr($sub,$i,$l);
					my $diffcount = () = ($que ^ $s) =~ /[^\x00]/g;
					if ($diffcount < $m_r){
						$jun = $i;
						$m_r = $diffcount;
					}
				}
				if($jun != -1){
					$jun = $jun + $pos_t;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTS\n";
					print STDERR "$rds\n";
				}else{
					$num_f ++;
				}
			}else{
				$ins_direc = "R";
				my $sub_start = $pos_t-$ins_len +$len_t;
				$sub = substr($genomes{$chr_t},$sub_start-1,$ins_len);
				my $jun = -1;
				my $m_r = 5;
				for my $i (0..(length ($sub)- length ($que))){
					my $s = substr($sub,$i,$l);
					my $diffcount = () = ($que_r ^ $s) =~ /[^\x00]/g;
					if ($diffcount < $m_r){
						$jun = $i;
						$m_r = $diffcount;
					}
				}
				if($jun != -1){
					$jun = $jun + $sub_start;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTS\n";
					print STDERR "$rds\n";	
			}else{
					$num_f++;
				}
			}
		}else{
			$num_f ++;
		}
	}
}

sub te_end{
	my %cors = @_;
	if ($cors{$te}{cig} =~ /E:(\d+)/){
		my $l = $1;
		my $end_pos = $cors{$te}{pos} + length($cors{$te}{seq}) - $l;
		if ($end_pos >= $tnt_len-10 or 	($end_pos <= ($ltr+5) and $end_pos >= ($ltr - 10))){
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
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_len);
				my $jun = -1;
				my $m_r = 2 ;
				for my $i (0..(length ($sub)- length ($que))){
					my $s = substr($sub,$i,$l);
					my $diffcount = () = ($que_r ^ $s) =~ /[^\x00]/g;
					if ($diffcount <= $m_r){
						$jun = $i;
						$m_r = $diffcount;
					}
				}
				if($jun != -1){
					$jun = $jun + $pos_t;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTE\n";
					print STDERR "$rds\n";
				}else{
					$num_f ++;
				}
			}else{
				$ins_direc = "S";
				my $sub_start = $pos_t+$len_t-$ins_len;
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_len);
				my $jun = -1;
				my $m_r = 2 ;
				for my $i (0..(length ($sub)- length ($que))){
					my $s = substr($sub,$i,$l);
					my $diffcount = () = ($que ^ $s) =~ /[^\x00]/g;
					if ($diffcount <= $m_r){
						$jun = $i;
						$m_r = $diffcount;
					}
				}
				if($jun != -1){
					$jun = $jun + $pos_t;
					print "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\tTE\n";
					print STDERR "$rds\n";
				}else{
					$num_f ++;
				}
			}
		}
	}
}

sub ge_start{
	my %cors = @_;
	if($cors{tar}{cig} =~ /S:(\d+)/){
		my $l = $1;
		my $que = substr($cors{tar}{seq},0,$l);
		( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
		$que_r = reverse $que_r; 
		my $len = (length $que ) ;
		my $sub_h = "AAAAA".substr($genomes{$te},0,$len+5);
		my $sub_t = substr($genomes{$te},-($len+5))."AAAAA";
		my $ma = match($que,$que_r,$sub_h,$sub_t);
		if ($ma != 0){
			my $ma_direc = $ma * $cors{tar}{direc};
			if (abs($ma_direc) == 3){
				my $ins_direc = "S";
				my $jun  = $cors{tar}{pos};
				print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGS\n";
				print STDERR "$rds\n";
				$num_eff++;
			}elsif(abs ( $ma_direc )== 2 ) {
				my $ins_direc = "R";
				my $jun = $cors{tar}{pos};
				print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGS\n";
				print STDERR "$rds\n";
				$num_eff++;
			}else{
				$num_f ++;
			}
		}else{
			$num_f ++;
		}
	}
}

sub ge_end{
	my %cors = @_;
	if($cors{tar}{cig} =~ /E:(\d+)/){
				my $l = $1;
				my $que = substr($cors{tar}{seq},-$l);
				( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
				$que_r = reverse $que_r; 
				my $len = (length $que ) ;
				my $sub_h = "AAAAA".substr($genomes{$te},0,$len+5);
				my $sub_t = substr($genomes{$te},-($len+5))."AAAAA";
				my $ma = match($que,$que_r,$sub_h,$sub_t);
				if ($ma != 0){
					my $ma_direc = $ma * $cors{tar}{direc};
					if (abs($ma_direc) == 3){
						my $ins_direc = "R";
						my $jun = $cors{tar}{pos} + length($cors{tar}{seq})-$l;
						print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGE\n";
						print STDERR "$rds\n";
						$num_eff++;
					}elsif ( abs($ma_direc) == 2 ){
						my $ins_direc = "S";
						my $jun = $cors{tar}{pos} + length($cors{tar}{seq})-$l;
						print "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGE\n";
						print STDERR "$rds\n";
						$num_eff++;
					}else{
						$num_f++;
					}
				}else{
					$num_f++;
				}
	}else{
			$num++;
	}
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
				#print "$diffcount\t$len_q\n";
				my $direc = ($j==0?1:-1);
				$relas{$diffcount} = $direc*$k;
			}
		}
	}
	foreach (sort {$a<=>$b} keys %relas){
		if ( $_ <= 3){
			my $v = $relas{$_};
			return "$v";   #### return value SH  ST RH  RT  
		}else{
			return "0";
		}
		last;
	}
}

