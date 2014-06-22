#!/usr/bin/perl
use warnings; use strict;
my($sample_name,$reads1,$reads2,$phred,$len,$ad_1,$ad_2,$file_sta) = @ARGV;
die "$0 <sample> <'reads1_files'> <'reads2_files'> <phred> <lenth> <5'adaptor:Uni Index> <3'adaptor> <full_path_stafile>\n" if(@ARGV != 8);

####read adapter###=

my %adaptors;
while (<DATA>){
	chomp;
	my ($id,$seq) = split /\t/,$_;
	$adaptors{$id} = $seq;
}
my $uni   = $adaptors{$ad_1};
$uni =~ tr/ATCGatcg/TAGCtagc/;
$uni = reverse $uni;

my $index = $adaptors{$ad_2};
$index = "A".$index;

print "Adaptor:\n$uni\n$index\n";
######## making your main result directory!!! ########
print "\n############making dir: $sample_name\n";
mkdir $sample_name;

########unzipping the original gz files########
print "\n############unzipping fastq files:\n\t$reads1\n\t$reads2\n.....\n";
my $fq_1 = $sample_name."_R1";
my $fq_2 = $sample_name."_R2";
system ("printf '$sample_name\tRaw:' >> $file_sta") == 0 or die $!;
system ("head -n 874749056 $reads1 |count_line_num.pl $file_sta >$sample_name/$fq_1") == 0  or die $!;
system ("head -n 874749056 $reads2 |count_line_num.pl $file_sta >$sample_name/$fq_2") == 0  or die $!;

####### change to your working dir############

chdir "$sample_name" or die $!;

#######    dynamic triming_step_1   ###################
my $dyn_dir = "Dynamic_dir_step_1";
mkdir $dyn_dir;
print "\n############dynamic triming reads1 $fq_1\n";
system ("DynamicTrim $fq_1 -h $phred -d $dyn_dir") == 0 or die $!;
print "###########dynamic triming reads2 $fq_2\n";
system ("DynamicTrim $fq_2 -h $phred -d $dyn_dir") == 0 or die $!;


my $fq_1_trimmed = "$dyn_dir/$fq_1.trimmed";
my $fq_2_trimmed = "$dyn_dir/$fq_2.trimmed";

system ("rm -f $fq_1") == 0 or die $!;
system ("rm -f $fq_2") == 0 or die $!;
############### Length_sorting from dynamic trim ###########
my $len_dir = "Len_sort_step_2";
mkdir $len_dir;
print "\n###########lenth sorting step 2############\n";
system ("LengthSort $fq_1_trimmed $fq_2_trimmed -l $len -d $len_dir") == 0 or die $!;
system ("rm -f $fq_1_trimmed $fq_2_trimmed") == 0 or die $!;


#######    cut adaptor   !!!########
my $cut_dir = "Cutadt_dir_step_3";
mkdir $cut_dir;
chdir $cut_dir;
my $fq_1_trimmed_cut = "$fq_1.trimmed.cut";
my $fq_2_trimmed_cut = "$fq_2.trimmed.cut";

print "\n#############cuting adapter file 1 $fq_1_trimmed\n";
my $fq_1_loc = "../$len_dir/$fq_1.trimmed.paired1";
my $fq_1_trimmed_cut_sta = "$fq_1_trimmed_cut.sta";
system ("printf 'Trim:' >>$file_sta") == 0 or die $!;
system ("cat $fq_1_loc | count_line_num.pl $file_sta | cutadapt -a $index  -f fastq - -o $fq_1_trimmed_cut > $fq_1_trimmed_cut_sta") == 0 or die $!;

print "\n#############cuting adapter file 2 $fq_2_trimmed\n";
my $fq_2_loc = "../$len_dir/$fq_1.trimmed.paired2";
my $fq_2_trimmed_cut_sta = "$fq_2_trimmed_cut.sta";
system ("cat $fq_2_loc | count_line_num.pl $file_sta | cutadapt -a $uni  -f fastq -  -o $fq_2_trimmed_cut > $fq_2_trimmed_cut_sta") == 0 or die $!;
chdir ".."; # change to your main result dir

system ("rm -f $len_dir/*.paired*") == 0 or die $!;
############## Length sort #########
my $len_dir_4 = "Len_sort_step_4";
mkdir "$len_dir_4";
print "\n#############length sorting\n";
system ("LengthSort $cut_dir/$fq_1_trimmed_cut $cut_dir/$fq_2_trimmed_cut -l $len -d $len_dir_4") == 0 or die $! ;
 
system ("printf 'Adaptor:' >>$file_sta") == 0 or die $!;
system ("cat $len_dir_4/$fq_1_trimmed_cut.paired1 | count_line_num.pl $file_sta NO") == 0 or die $!;
system ("cat $len_dir_4/$fq_1_trimmed_cut.paired2 | count_line_num.pl $file_sta NO") == 0 or die $!;
system ("echo  >>$file_sta") == 0 or die $!;

system ("rm -f $cut_dir/$fq_1_trimmed_cut") == 0 or die $!;
system ("rm -f $cut_dir/$fq_2_trimmed_cut") == 0 or die $!;
##########  compress file  ######
__DATA__
Uni	AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Index	GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
Index_1	GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
Index_2	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
Index_3	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
Index_4	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
Index_5	GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
Index_6	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
Index_7	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
Index_8	GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
Index_9	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
Index_10	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
Index_11	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
Index_12	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
Index_13	GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG
Index_14	GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG
Index_15	GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG
Index_16	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG
Index_18	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG
Index_19	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG
