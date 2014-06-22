#!/usr/bin/perl
use warnings; use strict;

my ($sam_file) = @ARGV;

sub read_sam{
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

my @aligns = read_sam($sam_file);
foreach my $reads (@aligns){
	print "isdfsd$reads\n";
}
