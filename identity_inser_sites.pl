#!/usr/bin/perl
use warnings; use strict;
use Seq;
use Getopt::Std;
use Bio::AlignIO;

########## get parameters #################
my %opts;
my $usage = "$0
	-h : help message
	-s : informative sam file
	-g : genome file
	-l : library insertion size
	-n : the id of te seq
	-a : <default : 0>, allow how many bases can be lost when transpositon,
	-p : 'directory/prefix' of your output files, relative to your working direcroty;
		 this script generate two output files : '\$WD/directory/prefix.ins.loc.lst' and '\$WD/directory/prefix.supported.reads.sam
	-r : tnt realian sam 
	-d : debug mode on 
";

die "$usage\n" if (@ARGV==0);
getopts("hs:g:l:n:a:p:r:d",\%opts);
if ($opts{h}){ print "$usage"; exit};

my ($sam_file,$genome_file,$ins_size,$te) = ($opts{s},$opts{g},$opts{l},$opts{n});
my $lost = $opts{a}?$opts{a}:0;
my $sam_te = $opts{r};
my $project = $opts{p};




############# put genome seq in hash  #########
my %genomes = Seq::seq_hash($genome_file);
my $tnt_len = length ($genomes{$te});  # transposon length



############# files used to save results ##
#open OUT, ">/dev/stdout" or die $!;
open OUT ,">${project}.$te.ins.loc.lst" or die $!;
open SUPP,">${project}.$te.support.reads.sam" or die $!;
open NSUPP,">${project}.$te.unsupport.reads.sam" or die $!;

############ parsing the te realn sam file ##### because TE have LTRs at both end

my %guanxi = te_aln($sam_te); 



###########################################################
################ the mainbody of code######################
###########################################################
###########################################################
my @aligns = scan_sam ($sam_file);

=head
foreach my $it ( @aligns){
	print "@$it\n" if $$it[0] =~ (/ACB053:89:D260YACXX:7:1107:12510:30746/);
}
=cut 

my $rds;
foreach my $grp (@aligns){  # parse each group of reads
	my %cors;
	my @hits = @$grp;
	$rds = join "\n",@hits;
	
	my %te_ha = read_2_ha($hits[0] );
	$cors{$te} = \%te_ha;
	my %chr_ha = read_2_ha($hits[1] );
	$cors{tar} = \%chr_ha;

	
	#print "THERE :$rds\n";
	my $tar_cig = $cors{tar}{cig};
	my $te_cig = $cors{$te}{cig}; 
	
	my $boo_ns = 0;
	if ( $tar_cig =~ /M/ and $te_cig =~ /M/){
		$boo_ns = 1 if( cross(%cors));
	}
	if($te_cig =~ /S/){
		$boo_ns = 1 if(te_start(%cors));
	}
	if($te_cig =~ /E/){
		$boo_ns = 1 if(te_end(%cors));
	}
	if($tar_cig =~ /S/){
		$boo_ns = 1 if (ge_start(%cors));
	}
	if($tar_cig =~ /E/){
		$boo_ns = 1 if (ge_end(%cors));
	}
	print NSUPP "$rds\n"  unless($boo_ns);
}
#################### the end of mainbody of code ######################




#######    sub functions  #######
sub te_aln{
	my $sam_te = shift @_;	
	
	my %guanxi;
	my %te_rcd;
	open my $fh, "samtools view -S -X $sam_te|"  or die $!;
	my $l_seq;
	my $l_as;
	while (<$fh>){
		chomp;
		my @ar  = (split /\t/,$_,12);
		$ar[5] =~ s/H/S/;
		my ($id,$flag,$chr,$pos,$cig,$seq,$tags) = @ar[0,1,2,3,5,9,11];
	

		my $direc = ($flag =~ /r/)?-1:1;
		my ($as) = $tags =~ /AS:i:(\d+)/;  

		## firstly, save the full infor to $l_seq and $l_as	
		unless($te_rcd{$id}){
			$guanxi{$id} = [];
			if ($direc == -1){
				$l_seq = Seq::rev_com($seq);
			}else{
				$l_seq = $seq;
			}
			$l_as = $as;
			$te_rcd{$id} = 1;
		}
		next if ($chr !~ $te);
		if($direc == -1){
			$seq = Seq::rev_com($l_seq);
		}else{
			$seq = $l_seq;
		}
		$ar[9] = "$seq";
		my $p = join "\t", @ar;
	

		if ($cig =~ /^\d+M$/){
			push @{$guanxi{$id}}, $p;
		}elsif($cig =~ /^\d+M(\d+)S$/){
			if (($tnt_len - (length($l_seq) - $1) +1 - $lost - 10) <= $pos){
				push @{$guanxi{$id}} ,$p;
			}
		}elsif($cig =~ /^\d+S\d+M$/){
			if ( $pos <= $lost	+ 10 ){
				push @{$guanxi{$id}} , $p;
			}
		}
	}
	return %guanxi;
}
	
	

sub scan_sam{    # put the pair reads in to one element of one array @re
	my $file = shift @_;
	open my $fh , $file or die $!;
	
	my @re;   
	while (<$fh>){
		chomp;
		if (/^@/){
			print SUPP "$_\n";
			next;
		}
		my($id,$flag,$chr) = (split /\t/,$_)[0,1,2];
		next if ( $chr =~ /$te/);

		my $r = ($flag =~ /1/)?1:2;
		$_ =~ s/$id\t/$id:$r\t/;
		my $r_a_t = ($r =~ /2/)?1:2;
		
		if( defined $guanxi{"$id:$r_a_t"}){
			foreach my $te_aln (@{$guanxi{"$id:$r_a_t"}}){
				#print "$te_aln\n$_\n" if ($id eq "SRR556175.40406632");
				push @re, [$te_aln,$_];
			}
		}
	}
	return @re;
}	

sub read_2_ha{     # this subroutine used to
	my %cors;
	my $hit =shift  @_;
	my($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$hit)[0,1,2,3,4,5,6,7,9];
	my $rc = (($flag =~ /r/)? -1:1);
	my $cs = cigar($cig);
	($id,my $r) = $id =~ /(.+)\:(\d)$/;
	$cors{cig} = $cs;
	$cors{direc} = $rc;
	$cors{id} = $id;
	$cors{pos} = $pos;
	$cors{seq} = $seq;
	$cors{chr} = $chr;
	$cors{mq} = $mq;
	#print "$cs\t$rc\t$id\n" if ($id =~ /SRR556175\.40406632/);
	return %cors;
}


sub cigar {
	my $cig = shift @_;
	return "Z" if ($cig =~ /\*/);
	my $len;		        #  reads length
	my $cs;                      #  simplified cigars 
	my ($m_l);					# match length
		
	while ($cig =~ /(\d+)([MSIH])/g){
		my $n = $1;
		my $c = $2;
		$len += $n;
		$m_l += $1 if ($c eq "M");
	}
	if ($m_l/$len > 95/100  ){     #  matcha and mismatch > 0.95 ; it will be considered as totally matched reads
		$cs = "M";
	}else{
		if ( $cig =~ /^(\d+)[SH]/ and $1 >= 5){
			$cs .= "S:$1";
		}
		if ( $cig =~ /(\d+)[SH]$/ and $1 >= 5){
			$cs .= "E:$1";
		}
	}
	$cs = ($cs?$cs:"Z");
	return $cs;
}

sub cross {
	my %cors = @_;
	my $len_tar = length ($cors{tar}{seq});
	my $len_te  = length ($cors{$te}{seq});
	my $len_inter = int($ins_size/2);
	if ( $cors{$te}{direc} == 1 and $cors{$te}{pos} > ($tnt_len - $ins_size  - $lost)){
		
		#   ---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
		#                                                               s----->    <--------
		#
		#   ----------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------
		#                 -------->       <------s
		my ($ins_direc,$jun);
		if ($cors{tar}{direc} == 1){
			$ins_direc = "R";
			$jun = $cors{tar}{pos} + $len_tar + $len_inter;  # assume the ins site is the end of match of reads at genome
		}else{
			$ins_direc = "S";
			$jun = $cors{tar}{pos} - $len_inter;                            # assume the ins site at the start of match of read at genome
		}
		my $te_jun = $cors{$te}{pos} + $len_te;
		print OUT "$cors{$te}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tCE:$te_jun\t$cors{tar}{mq}\n";
		print SUPP "$rds\n";
		return 1;
	}elsif ( $cors{$te}{direc} == -1  and $cors{$te}{pos} < ($ins_size  - $len_te + $lost )){
		
		#   ---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
		#                   -------->       s<--------
		#
		#   ----------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------
		#                                                               ------>s	<-----
       	my $len_tar = length ($cors{tar}{seq});	
		my ($ins_direc,$jun);
		if ( $cors{tar}{direc} == 1 ){
			$ins_direc = "S";
			$jun = $cors{tar}{pos} + $len_tar + $len_inter;
		}else{
			$ins_direc = "R";
			$jun = $cors{tar}{pos} - $len_inter;
		}
		my $te_jun = $cors{$te}{pos};
		print OUT "$cors{$te}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tCS:$te_jun\t$cors{tar}{mq}\n";
		print SUPP "$rds\n";
		return 1;
	}else{
		print STDERR  "error:cross fail:$cors{$te}{id}\tte_direc:$cors{$te}{direc}\tte_pos:$cors{$te}{pos}\ttar_direc:$cors{tar}{direc}\n" if ($opts{d});
	}
}

sub te_start{    # postion limitor in guanxi
	my %cors = @_;
	if($cors{$te}{cig} =~ /S:(\d+)/ and $cors{$te}{direc} == -1 ){
		
		# ---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
		#             ------->     <-------
	
		#
		# --------------------------- <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------
		#                                                                 ------->   <-----
		
		my $l = $1;
		
		my $que = uc(substr($cors{$te}{seq},0,$l));
		( my $que_r = $que) =~ tr/ATCG/TAGC/;
		$que_r = reverse $que_r;
		my $ins_direc;
		my $sub;
		my $chr_t = $cors{tar}{chr};
		my $pos_t = $cors{tar}{pos};
		my $len_t = length $cors{tar}{seq};
		my $jun_te = $cors{$te}{pos};
		if($cors{tar}{direc} == 1){
			$ins_direc = "S";
			$sub = substr($genomes{$chr_t},$pos_t-1,$ins_size);
			
			my ($diff,@juns) = mat($que,$sub);
			
			if($diff/$l < 0.05){
				#@juns = map { $_ + $pos_t+ $l -1 } @juns;
				#my $jun = join ":", @juns;
				my $jun = $juns[0] + $pos_t+ $l -1;
				my $mark = @juns >1 ?"ts":"TS";
				print OUT "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\t$mark:$jun_te\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}
		}else{
			$ins_direc = "R";
			my $sub_start = $pos_t-$ins_size +$len_t;
			$sub = substr($genomes{$chr_t},$sub_start-1,$ins_size);
			
			my($diff,@juns) = mat($que_r,$sub);
			if($diff/$l < 0.05){
				#@juns = map{ $_ + $sub_start} @juns ;
				#my $jun = join ":",@juns;
				my $jun = $juns[0] + $sub_start;
				my $mark = @juns >1 ?"ts":"TS";
				print OUT "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\t$mark:$jun_te\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}	
		}
	
	}else{
		print  STDERR "error: te_start_direc:$cors{$te}{id}\t$cors{$te}{cig}\t$cors{$te}{direc}\n" if ($opts{d});
	}
}

sub te_end{
	my %cors = @_;
	if ($cors{$te}{cig} =~ /E:(\d+)/ and $cors{$te}{direc} == 1){

		#---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-------------------------------
		#                                                      -------->     <-------
		#
		#----------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-------------------------------
		#                ------>  <--------      
		my $l = $1;
		my $end_pos = $cors{$te}{pos} + length($cors{$te}{seq}) - $l -1;
			
		my $que = uc(substr($cors{$te}{seq},-$l));
		( my $que_r = $que) =~ tr/ATCG/TAGC/;
		$que_r = reverse $que_r;
		my $ins_direc;
		my $sub;
		my $chr_t = $cors{tar}{chr};
		my $pos_t = $cors{tar}{pos};
		my $len_t = length $cors{tar}{seq};
		if($cors{tar}{direc} == 1){
			$ins_direc = "R";
			$sub = substr($genomes{$chr_t},$pos_t-1,$ins_size);
			
			my ($diff,@juns) = mat($que_r,$sub);
			
			if( $diff/$l <0.05 ){
				#@juns = map {$_ + $pos_t + $l -1} @juns;
				#my $jun = join ":",@juns;
				my $jun = $juns[0] + $pos_t + $l -1;
				my $mark = @juns > 1 ? "te":"TE";
				print OUT "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\t$mark:$end_pos\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}
		}else{
			$ins_direc = "S";
			my $sub_start = $pos_t+$len_t-$ins_size;
			$sub = substr($genomes{$chr_t},$sub_start-1,$ins_size);
			
			my ($diff,@juns) = mat($que,$sub);
			if($diff/$l < 0.05 ){
				#@juns = map {$_ + $sub_start} @juns;
				#my $jun = join ":",@juns; 
				my $jun = $juns[0] + $sub_start;
				my $mark = @juns > 1 ? "te":"TE";
				print OUT "$cors{$te}{id}\t$ins_direc\t$chr_t\t$jun\t$mark:$end_pos\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}
		}
	}else{
		print "error  te_end_direc:$cors{$te}{id}\t$cors{$te}{cig}\t$cors{$te}{direc} \n"if ($opts{d});
	}
}



sub ge_start{
	my %cors = @_;
	if($cors{tar}{cig} =~ /S:(\d+)/ ){
		
		my $l = $1;
		
		my $que = uc(substr($cors{tar}{seq},0,$l));
		my $ry = uc(substr($cors{tar}{seq},$l));
		( my $que_r = $que) =~ tr/ATCG/TAGC/;
		( my $ry_r = $ry) =~ tr/ATCG/TAGC/;
		$que_r = reverse $que_r;	
		$ry_r  = reverse $ry_r;

		my $sub_h = "NNNNNNNNNN".substr($genomes{$te},0,$l+$lost+10);
		my $sub_t = substr($genomes{$te},-($l+$lost+10))."NNNNNNNNNN";
		
		my %ma_hash = match($que,$que_r,$sub_h,$sub_t);	
		if  (%ma_hash){
			if ($ma_hash{3}  and $cors{$te}{direc} == 1 and $cors{$te}{pos} > ($tnt_len - $ins_size - $lost - 10 )  ){
				
		#  --------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-------------------------------------------------
		#                                                         -------->        <---------
		#
				my $adj = $ma_hash{3}[-1];
				my $seq_te = substr($sub_t,$adj+$l);
				my $r = 0;
				for my $i (0..length($seq_te)){
					my $s = substr($seq_te,$i,1);
					my $g = substr($ry,$i,1);
					if ($s eq $g){
						$r = $i +1; 
					}else{
						last;
					}
				}
				my $rr = length($seq_te) - 10 - $r;
				#my $mat = substr($sub_t,$adj,$l);
				#print OUT "RAW:$mat\t$seq_te\t$que\t$ry\n";	
				my $ins_direc = "S";
				my $jun  = $cors{tar}{pos} + $r;
				my $te_jun = $tnt_len - $rr;
				#my $te_jun = $tnt_len -$lost + $adj + $r;
				print OUT "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGE:$te_jun\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}elsif($ma_hash{-2}  and $cors{$te}{direc} == -1 and $cors{$te}{pos} < ($ins_size - length($cors{$te}{seq}) + $lost + 10)) {
				
		#
		#  ---------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-------------------------------------------------
		#                                                         -------->        <---------  
				my $adj =$ma_hash{-2}[0];
				my $seq_te = substr($sub_h,0,$adj);
			
				#my $mat = substr($sub_h,$adj,$l);
				#print OUT "RAW:$mat\t$seq_te\t$que_r\t$ry_r\n";
				
				my $r = 0;
				for my $i (0..length($seq_te)){
					my $s = substr($seq_te,-($i+1));
					my $g = substr($ry_r,-($i+1));
					if($s eq $g){
						$r = $i+1;
					}else{
						last;
					}
				}
				my $rr = length($seq_te) - 10 - $r;
				my $ins_direc = "R";
				my $jun = $cors{tar}{pos} + $r ;
				my $te_jun = $rr +1;	
				#my $te_jun = $adj + $rr -10 +1;
				print OUT "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGS:$te_jun\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}
		}else{
			print "error ge_start_mism:$cors{tar}{id}:$que\t$que_r\n"if ($opts{d});
		}
	}
}

sub ge_end{
	my %cors = @_;
	if($cors{tar}{cig} =~ /E:(\d+)/){
	
		my $l = $1;
		my $tar_l = length($cors{tar}{seq});
				
		my $que = uc(substr($cors{tar}{seq},-$l));
		my $ry = uc(substr($cors{tar}{seq},0,$tar_l-$l));
		( my $que_r = $que) =~ tr/ATCG/TAGC/;
		( my $ry_r = $ry)  =~ tr/ATCG/TAGC/;
		$que_r = reverse $que_r; 
		$ry_r = reverse $ry_r;

		my $sub_h = "NNNNNNNNNN".substr($genomes{$te},0,$l+$lost+10);
		my $sub_t =  substr($genomes{$te},-($l+$lost+10))."NNNNNNNNNN";

		my %ma_hash = match($que,$que_r,$sub_h,$sub_t);

		if (%ma_hash){
			
		#
		#  ---------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-------------------------------------------------
		#                                     -------->        <----------  
			

			if ( $ma_hash{-3} and $cors{$te}{direc} == 1 and $cors{$te}{pos} > ($tnt_len - $ins_size - $lost -10 )){
				my $adj = $ma_hash{-3}[-1];
				my $seq_te = substr($sub_t,$adj+$l);
				my $r = 0;
				for  my $i (0..length($seq_te)){
					my $s = substr($seq_te,$i,1);
					my $g = substr($ry_r,$i,1);
					if ($s eq $g){
						$r = $i +1;
					}else{
						last;
					}
				}
				#my $mat = substr($sub_t,$adj,$l);
				#print OUT "RAW:$mat\t$seq_te\t$que_r\t$ry_r\n";
				my $rr = length($seq_te) - 10 - $r;
				my $ins_direc = "R";
				my $jun = $cors{tar}{pos} + (length($cors{tar}{seq})-$l) - $r -1 ;
				my $te_jun = $tnt_len - $rr; 
				#my $te_jun = $tnt_len -$lost +$adj + $rr;

				print OUT "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGE:$te_jun\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}elsif ( $ma_hash{2} and $cors{$te}{direc} == -1 and $cors{$te}{pos} < ($ins_size - length($cors{$te}{seq}) + $lost + 10)){
				my $adj = $ma_hash{2}[0];

		#  --------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-------------------------------------------------
		#                                    --------->       <---------
		#
				
				my $seq_te = substr($sub_h,0,$adj);
				
				#my $mat = substr($sub_h,$adj,$l);
				#print OUT "RAW:$mat\t$seq_te\t$que\t$ry\n";
	
				my $r = 0 ;
				for my $i (0..length($seq_te)){
					my $s = substr($seq_te,-($i+1));
					my $g = substr($ry,-($i+1));
					if($s eq $g){
						$r = $i+1;
					}else{
						last;
					}
				}
				my $rr = length($seq_te) - $r -10;

				my $ins_direc = "S";
				my $jun = $cors{tar}{pos} + (length($cors{tar}{seq}) - $l) - $r -1;
				my $te_jun = $rr +1;
				#my $te_jun = $adj + $rr -10 +1;

				print OUT "$cors{tar}{id}\t$ins_direc\t$cors{tar}{chr}\t$jun\tGS:$te_jun\t$cors{tar}{mq}\n";
				print SUPP "$rds\n";
				return 1;
			}
		}else{
			print "error: ge_end_mism:$cors{tar}{id}\t$que\tsub:$sub_t\t$sub_h\n"if ($opts{d});
		}
	}
}

sub mat { 
	my ($que,$sub) = @_;
	$que = uc($que);
	$sub = uc($sub);
	my $q_l = length $que;
	my $s_l = length $sub;
	
	my %ref_dif;
	for my $i (0..($s_l-$q_l)){
		my $tgt = substr($sub,$i,$q_l);
		my $diffcount = () = ($que ^ $tgt) =~ /[^\x00]/g;
		$ref_dif{$diffcount} = [] unless(exists $ref_dif{$diffcount});
		push @{$ref_dif{$diffcount}}, $i;
	}
	foreach my $k (sort {$a <=> $b} keys %ref_dif){
		my @pos = @{$ref_dif{$k}};
		return ($k, @pos);   #  0 based return position
		last;
	}
}


# return one hash with key if direction and value is start position of match on subject
sub match {
	my %relas;
	
	for my $j (0..1){
		for my $k (2..3){
			my $que = $_[$j];
			my $sub = $_[$k];
			my ($diff,@locs) = mat($que,$sub);
			my $ratio = $diff/(length $que);
			my $direc = ($j==0?1:-1);
			
			my $dir_pos  = $direc * $k;	
			$relas{$ratio}{$dir_pos} = [@locs]   # head of tail of sequence matched
		}
	}
	
	# exact the most similarity region and must be small than 0,05
	foreach (sort {$a<=>$b} keys %relas){	
		if ( $_ <= 0.05){
			my $v = $relas{$_};
			my %m_h = %$v;
			return (%m_h);   #### return value SH  ST RH  RT  
		}
		last;
	}
}

sub count_m {
	my $cig = shift @_;
	my $num;
	while($cig =~ /(\d+)M/g){
		$num += $1;
	}
	return $num;
}
