#!/usr/bin/perl
use warnings; use strict;

my ($ref,@rest) = @ARGV;

open REF, $ref or die $!;

my %ref;  # save each point in hash, value is the insert position
my %pri;  # key is the insert position, value is the support information
my @rec;  # save the order of insert position;
while(<REF>){
	chomp;
	my ($chr,$s,$e,$t) = split /\t/, $_;

	foreach my $i ( $s..$e){
		$ref{$chr}{$i} = "$chr:$s:$e";
	}
	$pri{"$chr:$s:$e"} = $t;
	push @rec,"$chr:$s:$e";
}

close(REF);

my %other;
for my $i (@rest){   # process other file
	open OT, $i or die $!;
	my %tem_ha;    # temperate hash for save information for individual file
	while(<OT>){
		chomp;	
		my ($c,$s,$e,$t) = split /\t/, $_;
		my $site;   # the position in the ref list
		for my $j ($s..$e){
			if (exists $ref{$c}{$j}){
				$site = $ref{$c}{$j};
				last;
			}
		}
		$tem_ha{$site} = "$t" if ($site);   # correspond the type with site in ref
	}
	$other{$i} = \%tem_ha;     # save in %other hash
}

my @list;
print (join "\t",("POS",$ref,@rest,"\n") );
for my $i (@rec){    #  $i is the insert position in ref file
	my $p = "$i\t$pri{$i}";     # the support record in ref file
	for my $j (@rest){      # extract support record in other files
		my %ha = %{$other{$j}};   # dereference of hash
		$p .= $ha{$i}? "\t$ha{$i}":"\tNA";		# print
	}
	print "$p\n";
}

