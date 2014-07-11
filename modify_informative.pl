#!/usr/bin/perl
use warnings; use strict;
my ($diff_file) = @ARGV;
use Seq;

my %reads;
open my $fh, "samtools view -h -S -X $diff_file|" or die $!;

while (<$fh>){
	my @pairs;
	chomp;
	if (/^@/){
		print "$_\n";
		next;
	}
	my ($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$seq) = (split /\t/,$_)[0,1,2,3,4,5,6,7,9];
	if (! keys %reads or exists $reads{$id}){
		push @{$reads{$id}},$_; 
	}else{
		my ($g) = values %reads;
		my @hits = @$g;
		my %seq_ha;
		my %qua_ha;

		################ main body ###########
		foreach my $hit (@hits){
			my ($id,$flag,$chr,$pos,$mq,$cig,$nchr,$npos,$dis,$seq,$qua ) = (split /\t/,$hit)[0,1,2,3,4,5,6,7,8,9,10];
			my $r = ($flag =~ /1/?1:2);
			if ($cig !~ /H/){
				my ($s,$q);
				if($flag =~ /r/){
					$s = Seq::rev_com($seq);
					$q = reverse($qua);
				}else{
					$s = $seq;
					$q = $qua;
				}
				$seq_ha{$r} = $s;
				$qua_ha{$r} = $q;
				print "$hit\n";
			}else{
				my $s;
				my $q;
				if ($flag =~ /r/){
					$s = Seq::rev_com($seq_ha{$r});
					$q = reverse($qua_ha{$r});
				}else{
					$s = $seq_ha{$r};
					$q = $qua_ha{$r};
				}
				$cig =~ s/H/S/g;
				print "$id\t$flag\t$chr\t$pos\t$mq\t$cig\t$nchr\t$npos\t$dis\t$s\t$q\n";
			}
		}




		################## END ####################
		undef(%reads);
		push @{$reads{$id}}, $_;

	}
}


