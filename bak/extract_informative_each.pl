#!/usr/bin/perl
use warnings; use strict;


my @words = @ARGV;

my @fhs;
foreach (@words){
	open my $f,">$_.sam" or die $!;
	push @fhs,$f;
}

my %reads;
while (<STDIN>){
	chomp;
	next if (/^@/);
	my ($title,$chr,$cig,$rnext) = (split /\t/,$_)[0,2,5,6];
	if (! keys(%reads )  or exists ($reads{$title})){
		$reads{$title} .= "############$_";
	}else{
		my @arr = values %reads;
		my @va = split '############',$arr[0];
		shift @va;
		my $pt = join "\n",@va;
		LOOP: foreach my $it (@va){
			my ($title,$flag,$chr,$cig,$rnext) = (split /\t/,$it)[0,1,2,5,6];
			for (my $i = 0; $i < @words;$i++){
				if(($rnext =~ /$words[$i]/ and $chr !~ /$words[$i]/) or ($rnext ne "=" and $chr =~ /($words[$i])/)){
					my $fh = $fhs[$i];
					print $fh  "$pt\n";
					last LOOP;
				}
			}
		}
		undef(%reads);
		$reads{$title} = "############$_";
	}
}

