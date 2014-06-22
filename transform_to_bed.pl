#!/usr/bin/perl
use warnings; use strict;

my($ins_file,$window ) = @ARGV;

open INS,$ins_file or die $!;
chomp (my @sites = <INS>);
for (my $i = 0;;){
	my @clu;
	my $step ;
	push @clu,$sites[$i];

	###############  detect a strint of pos at genome #################
	for (my $j =1;;$j++){
		$step = $j;
		my $pre = $i + $j -1;
		my $nex = $i + $j;
		last if ($nex >= @sites);	
		my ($id_p,$dir_p,$chr_p,$pos_p,$ty_p) = split /\t/,$sites[$pre];
		my ($id_n,$dir_n,$chr_n,$pos_n,$ty_n) = split /\t/,$sites[$nex];
		if (($dir_p eq $dir_n ) and  ($chr_p eq $chr_n ) and  (($pos_n - $pos_p ) < $window)){
			push @clu,$sites[$nex];
		}else{
			last;
		}

	}

	############### if have more than 5 pos, It is a real insertion site ###
	if (@clu > 5){
		my %sit_s;
		my %sit_e;
		my %sit_r;
		my($id,$dir,$chr,$pos,$ty);
		foreach my $pos (@clu){
			($id,$dir,$chr,$pos,$ty) = split /\t/,$pos;
			if ($ty =~ /GS|TS/ ){
				$sit_s{$pos} ++;
			}elsif ($ty =~ /GE|TE/){
				$sit_e{$pos} ++;
			}
			$sit_r{$pos} ++;
		}
		my ($s) = (sort {$sit_s{$a}<=>$sit_s{$b}} keys %sit_s)[-1] ;
		my ($e) = (sort {$sit_e{$a}<=>$sit_e{$b}} keys %sit_e)[-1];
		my ($s_f,$e_f) = (sort {$sit_r{$a}<=>$sit_r{$b}} keys %sit_r)[0,-1];
		my $sc = scalar @clu;
		my $clu_pre = join "\n",@clu;
		#print "$i:$step\n$clu_pre\n";
		if ($s and $e){
			my ($ss,$ee ) = sort {$a<=>$b} ($s,$e);
			$ss --;
			print "$chr\t$ss\t$ee\t$sc:$dir:P\n";

		}else{
			my ($ss,$ee ) = sort {$a<=>$b}($s_f,$e_f);
			$ss --;
			print "$chr\t$ss\t$ee\t$sc:$dir:N\n";
		}
	}else{
		#print "$i:$step\n";
	}
	



	###############   end   #########################

	$i = $i + $step;
	last if ($i >= @sites);
}


			

