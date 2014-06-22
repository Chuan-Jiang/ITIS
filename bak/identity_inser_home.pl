#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Bio::Seq;
my ($sam_file,$genome_file,$ins_len,$te,$ltr) = @ARGV;
print "$0 <sam> <genome> <ins_len> <te_name> <ltr_len>\n" unless @ARGV;
############# genome seq in hash  #########
my %genomes = Seq::seq_hash($genome_file);
#print %simu_ge;
my $tnt_len = length ($genomes{$te});


###########################################
   ## Scan sam file ##

print STDERR "scanning sam file.....\n";
open SAM, $sam_file or die $!;
my %reads;
my $num;
my $num_eff;
my $num_f;
my $disc;
while(<SAM>){
	chomp;
	my($id) = (split /\t/,$_)[0];
	if(! keys %reads or exists ($reads{$id})){
		$reads{$id} .= "############$_";
	}else{
		my @arr = values %reads;
		my @va = split '############',$arr[0];
		shift @va;
		find(@va);
		undef(%reads);
		$reads{$id} = "############$_";
	}
}
print "$disc\n";
##################     parsing the paired reads
my $tnt;
my $tar;
my $rds;
sub find{
	$num ++;
	my @reads = @_;
	$rds = join "\n",@_; ## in order to print informative reads in one file 
	
	my %cors;
	my $cors_ref = \%cors;
	foreach my $it (@reads){
		my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$it)[0,1,2,3,4,5,6,7,9];	
		(my $r = $flag ) =~ s/\w+(\d)/$1/;  #  which reads?????
		if ($cig =~ /H/){
			$cors{$r}{hard} = $it;
			next;
		}			#  get rid of secondary alignment
		my $rc = (($flag =~ /r/)? -1:1);    # -1 antisense   and 1 means sense strand
		my $len = length $seq;		        #  reads length

		my $c;                      #  simplified cigar 
		my ($m);					# match length
		my $s = 0;					# match the longest soft clipped sequence
		while ($cig =~ /(\d+)([MS])/g){
			$m += $1 if ($2 eq "M");
			$s = $1 if ($2 eq "S" and $1 > $s);
		}
		if ($m/$len > 0.95 ){
			$c = "M";
		}elsif ( $cig =~ /^${s}S/ and $s >= 5){
			$c = "S:$s";
		}elsif ( $cig =~ /${s}S$/ and $s >= 5){
			$c = "E:$s";
		}else{
			$c = "Z";
		}
		
		$cors{$r}{cig} = $c;      # simplified cigar: M H T Z
		$cors{$r}{direc} = $rc;   # reads direction
		$cors{$r}{chr} = $chr;    # which chromosome
		$cors{$r}{id} = $id;	  # read id 
		$cors{$r}{pos} = $pos;  ### when reads if reverse complement the start loc of reads should plus the length of reads
		$cors{$r}{seq} = $seq;
		$cors{$r}{mq} = $mq;
	}

	if (($cors{1}{chr} =~ /$te/ )and  ($cors{2}{chr} ne $cors{1}{chr})){ 	
		$tnt = 1;
		$tar = 2;
	}elsif(($cors{2}{chr} =~ /$te/) and  ($cors{1}{chr} ne $cors{2}{chr})){ 
		$tnt = 2;
		$tar= 1;
	}elsif($cors{1}{chr} =~ /$te/ and $cors{2}{chr} =~ /$te/){ 
		$tnt = 3;
		$tar = 0;
	}elsif ($cors{1}{chr} !~ /$te/ and $cors{2}{chr} !~ /$te/ and ($cors{1}{chr} eq $cors{2}{chr})){
		$tnt = 0 ;
		$tar = 3;
	}
	
	if ($tnt > 0 and $tar >0){
		if ( $cors{$tar}{cig} eq "M" and $cors{$tnt}{cig} eq "M"){
			cross(%cors);
		}elsif($cors{$tnt}{cig} =~ /S/){
			te_start(%cors);
		}elsif($cors{$tnt}{cig} =~ /E/){
			te_end(%cors);
		}elsif($cors{$tar}{cig} =~ /S/){
			ge_start(%cors);
		}elsif($cors{$tar}{cig} =~ /E/){
			ge_end(%cors);
		}
		#	print "$tnt\t$tar\n";
	}elsif($tar == 3){

	}elsif($tnt == 3){

	}
	$tnt =0 ;
	$tar =0;

}


			






sub cross {
	my %cors = @_;
	if ( $cors{$tnt}{direc} == 1 and ( $cors{$tnt}{pos} > ($tnt_len - $ltr) or $cors{$tnt}{pos} < ($ltr - 100))){
		my ($ins_direc,$jun);
		if ($cors{$tar}{direc} == 1){
			$ins_direc = "R";
			$jun = $cors{$tar}{pos} + length ($cors{$tar}{seq});
		}else{
			$ins_direc = "S";
			$jun = $cors{$tar}{pos};
		}
		print "$cors{$tnt}{id}\t$ins_direc\t$cors{$tar}{chr}\t$jun\tCE\n";
		print STDERR "$rds\n";
		$num_eff ++;
	}elsif ( $cors{$tnt}{direc} == -1 and ($cors{$tnt}{pos} <  ($ltr-100)  or $cors{$tnt}{pos} > ($tnt_len - $ltr) )){
		my ($ins_direc,$jun);
		if ( $cors{$tar}{direc} == 1 ){
			$ins_direc = "S";
			$jun = $cors{$tar}{pos} + length($cors{$tar}{seq});
		}else{
			$ins_direc = "R";
			$jun = $cors{$tar}{pos};
		}
		print "$cors{$tnt}{id}\t$ins_direc\t$cors{$tar}{chr}\t$jun\tCS\n";
		print STDERR "$rds\n";
		$num_eff++;
	}else{
		$num_f ++;
	}
}

sub te_start{
	my %cors = @_;
	if($cors{$tnt}{cig} =~ /S:(\d+)/ ){
		my $l = $1;
		if ($cors{$tnt}{pos} <=2 or ($cors{$tnt}{pos} <= ($tnt_len - $ltr + 2) and $cors{$tnt}{pos} > ($tnt_len - $ltr))){
			my $que = substr($cors{$tnt}{seq},0,$l);
			( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
			$que_r = reverse $que_r;
			my $ins_direc;
			my $sub;
			my $chr_t = $cors{$tar}{chr};
			my $pos_t = $cors{$tar}{pos};
			my $len_t = length $cors{$tar}{seq};
			if($cors{$tar}{direc} == 1){
				$ins_direc = "S";
				$sub = substr($genomes{$chr_t},$pos_t-1,$ins_len);
				my $jun = -1;
				my $m_r = 2 ;
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
					print "$cors{$tnt}{id}\t$ins_direc\t$chr_t\t$jun\tTS\n";
					print STDERR "$rds\n";
				}else{
					$num_f ++;
				}
			}else{
				$ins_direc = "R";
				my $sub_start = $pos_t+$len_t-$ins_len;
				$sub = substr($genomes{$chr_t},$sub_start-1,$ins_len);
				my $jun = -1;
				my $m_r = 2;
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
					print "$cors{$tnt}{id}\t$ins_direc\t$chr_t\t$jun\tTS\n";
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
	if ($tnt > 0 and $cors{$tnt}{cig} =~ /E:(\d+)/){
		my $l = $1;
		my $end_pos = $cors{$tnt}{pos} + length($cors{$tnt}{seq}) - $l;
		if ($end_pos >= $tnt_len-10 or 	($end_pos <= ($ltr+5) and $end_pos >= ($ltr - 10))){
			my $que = substr($cors{$tnt}{seq},-$l);
			( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
			$que_r = reverse $que_r;
			my $ins_direc;
			my $sub;
			my $chr_t = $cors{$tar}{chr};
			my $pos_t = $cors{$tar}{pos};
			my $len_t = length $cors{$tar}{seq};
			if($cors{$tar}{direc} == 1){
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
					print "$cors{$tnt}{id}\t$ins_direc\t$chr_t\t$jun\tTE\n";
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
					print "$cors{$tnt}{id}\t$ins_direc\t$chr_t\t$jun\tTE\n";
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
	if($cors{$tar}{cig} =~ /S:(\d+)/){
		my $l = $1;
		my $que = substr($cors{$tar}{seq},0,$l);
		( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
		$que_r = reverse $que_r; 
		my $len = (length $que ) ;
		my $sub_h = "AAAAA".substr($genomes{$te},0,$len+5);
		my $sub_t = substr($genomes{$te},-($len+5))."AAAAA";
		my $ma = match($que,$que_r,$sub_h,$sub_t);
		if ($ma != 0){
			my $ma_direc = $ma * $cors{$tar}{direc};
			if (($ma_direc == -3 and $cors{$tnt}{direc} == 1)){
				my $ins_direc = "S";
				my $jun  = $cors{$tar}{pos};
				print "$cors{$tar}{id}\t$ins_direc\t$cors{$tar}{chr}\t$jun\tGS\n";
				print STDERR "$rds\n";
				$num_eff++;
			}elsif(( $ma_direc == 2 and $cors{$tnt}{direc} == -1)) {
				my $ins_direc = "R";
				my $jun = $cors{$tar}{pos};
				print "$cors{$tar}{id}\t$ins_direc\t$cors{$tar}{chr}\t$jun\tGS\n";
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
	if($tar >0 and $cors{$tar}{cig} =~ /E:(\d+)/){
				my $l = $1;
				my $que = substr($cors{$tar}{seq},-$l);
				( my $que_r = $que) =~ tr/ATCGatcg/TAGCtagc/;
				$que_r = reverse $que_r; 
				my $len = (length $que ) ;
				my $sub_h = "AAAAA".substr($genomes{$te},0,$len+5);
				my $sub_t = substr($genomes{$te},-($len+5))."AAAAA";
				my $ma = match($que,$que_r,$sub_h,$sub_t);
				if ($ma != 0){
					my $ma_direc = $ma * $cors{$tar}{direc};
					if (($ma_direc == -3 and $cors{$tnt}{direc} == 1)){
						my $ins_direc = "R";
						my $jun = $cors{$tar}{pos} + length($cors{$tar}{seq})-$l;
						print "$cors{$tar}{id}\t$ins_direc\t$cors{$tar}{chr}\t$jun\tGE\n";
						print STDERR "$rds\n";
						$num_eff++;
					}elsif ( $ma_direc == 2 and $cors{$tnt}{direc} == -1){
						my $ins_direc = "S";
						my $jun = $cors{$tar}{pos} + length($cors{$tar}{seq})-$l;
						print "$cors{$tar}{id}\t$ins_direc\t$cors{$tar}{chr}\t$jun\tGE\n";
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

