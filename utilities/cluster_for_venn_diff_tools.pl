#!/usr/bin/perl
use warnings; use strict;


my %rcd;
while(<>){
	chomp;
	my($chr,$s,$t) = (split /\t/, $_)[0,1,3];
	if($t =~ /^(\S+?),/){
		$t = $1;
	}

	$t = "ITIS" if ($t =~ /SR/);
	$rcd{$chr}{$s}{$t} = 1;
}

my %lst;
foreach my $chr(keys %rcd){
	my $st = $rcd{$chr};
	my %sub_ha = %$st;

	my $l_k=0;
	foreach my $k( sort {$a<=>$b} keys %sub_ha){   # $k contain the name of chromosome
		my $v = $sub_ha{$k};
		my @ts = keys (%$v);
		if ($k - $l_k  < 1000){
			foreach my $t (@ts){
				push @{$lst{$t} }, "$chr:$l_k";
			}
		}else{
			foreach my $t ( @ts){
				push @{$lst{$t} }, "$chr:$k";
			}
			$l_k = $k;
		}
	}
}

foreach my $tool ( keys %lst){
	my @locs = @{$lst{$tool}};
	my $num = @locs;
	print STDERR "#$tool\t$num\n";
	print "$tool\t".(join "\t", @locs) . "\n";
}


