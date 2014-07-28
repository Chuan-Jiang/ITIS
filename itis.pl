#!/usr/bin/perl
use warnings; use strict;
use FindBin;
use Getopt::Std;
use File::Basename;
use Time::localtime;

#########################  parameters #####################
my $usage = "USAGE:
	$0
	REQUIRED -g the genome sequence file in fasta format  
	REQUIRED -t the TE sequence file in fasta format
	REQUIRED -l the average length of fragments in library	
	REQUIRED -N the name of you project
	REQUIRED -1 the paired read file 1
	REQUIRED -2 the paired read file 2	
	REQUIRED -f gff file
			 
		-F <Y|N: default N> run scripts in 'FAST' mode; It won't align all reads to reference genome and caculate the bg depth 
			parameters specific to  '-F N'  :
			-B use your previous sorted and indexed bam file of all clean reads align to reference genome; on condition of '-F N'
			-R the minimum ratio of support fragments  and depth of 200 bp around the insertion site; default 0.2
			-D <3,200>, the depth range to filter candidate insertion site. 

		-q default: 1  the minimum average mapping quality of all supporting reads

		-e <Y|N: default N> if TE sequence have homolog in genome. using blast to hard mask repeat sequence is required
		
		-a <10> the allow number of base can be lost during transposon

		-b the total number of required supporting reads [3], minimum required supporting reads cover TE start[1], reads cover TE end[1]
			defualt: in the form of '3,1,1', seperate by comma 
		
		-c FORMAT:\\d,\\d,\\d; for  cpu number for 'BWA mem', 'samtools view'  and 'samtools sort'    defualt 8,2,2
		
		-w window size for cluster you support reads: DEFAULT :100bp
		
		-T use this specifed temperate directory or use the DEFAULT one :[project].[aStringOfNumbers]
		
		-m <Y|N: default F> Only print out all commands to STDERR
	         
		-h print STDERR this help message    

	
		eg: perl $0 -g genome.fa -t tnt1.fa -l 300 -n tnt1 -N test_run -1 reads.fq1 -2 reads.fq2 -f medicago.gff3 

		BWA samtools should in you PATH
	
	";

die "$usage\n" if (@ARGV == 0);
my %opt;
getopts("g:t:l:N:1:2:f:b:R:B:D:c:q:e:a:F:T:w:m:h",\%opt);

die "$usage\n" if ($opt{h});

########===============================

my $genome = $opt{g};
my $te_seq = $opt{t};
	my $te_base = basename $te_seq;
my $lib_len = $opt{l};
my $proj   = $opt{N};
my $rs1_ori  = $opt{1};
my $rs2_ori  = $opt{2};
my $gff = $opt{f};

my $min_reads  = $opt{b}?$opt{b}:"3,1,1";
my $ratio = $opt{R}?$opt{R}:"0.2";
my $bam = $opt{B}?$opt{B}:0;
my $depth_range= $opt{D}?$opt{D}:"3,200";
my $cpu    = $opt{c}?$opt{c}:"8,2,2";
	my($cpu_bwa,$cpu_view,$cpu_sort) = split /,/,$cpu;
my $exists = $opt{e}?$opt{e}:"N";
my $fast   = $opt{F}?$opt{F}:"N";
my $tmp_dir = $opt{T}? $opt{T} : "$proj.".time();
my $window = $opt{w}?$opt{w}:100;
my $only_cmd = $opt{m}?$opt{m}:"N";
my $cmd;
my $bindir  = "$FindBin::Bin";
my $map_q = $opt{q}?$opt{q}:1;
my $lost = $opt{a}?$opt{a}:10;
##########################################################
  

if(-e $tmp_dir){
	print "using dir: $tmp_dir\n";
}else{
	$cmd = "mkdir $tmp_dir";
	process_cmd($cmd);					# temperate folder contain intermidiate files
}

####################### prepare template #########

my $para_filter;
if($exists =~ /N/i){
	$para_filter = "";
	$cmd = "cat $genome $te_seq >$tmp_dir/$proj.ref_and_te.fa";
}else{
	$cmd = "perl $bindir/mask_te_homo_in_genome.pl -g $genome -t $te_seq -p $tmp_dir/te_homo_in_ref.lst -o $tmp_dir/$proj.ref_and_te.fa";
	$para_filter = "-l $tmp_dir/te_homo_in_ref.lst";
}


process_cmd($cmd);				# cat sequence together


$cmd = "cp $te_seq $tmp_dir/";
process_cmd($cmd);				# copy te sequence to tmp/


###################### Index you sequence file #######


$cmd = "bwa index $tmp_dir/$proj.ref_and_te.fa";
if ( -e "$tmp_dir/$proj.ref_and_te.fa.bwt"){
	print STDERR "Seems the Indexes for merged sequence exists. Skipped\n";
}else{
	process_cmd($cmd);				# index merged sequence
}


$cmd = "bwa index $tmp_dir/$te_base";
process_cmd($cmd);				# index te sequence


##### align original reads to reference genome ######
my $transformtobed_bam ;
if($fast =~ /N/i and $bam == 0){
	$cmd = "bwa mem -T 20 -t $cpu_bwa $tmp_dir/$proj.ref_and_te.fa $rs1_ori $rs2_ori | samtools view -@ $cpu_view -bS - | samtools sort -@ $cpu_sort - $tmp_dir/$proj.aln.bg.ref.sort";
	process_cmd($cmd);
	$transformtobed_bam = "-b $tmp_dir/$proj.aln.bg.ref.sort.bam";
	
	$cmd = "samtools index $tmp_dir/$proj.aln.bg.ref.sort.bam";
	process_cmd($cmd);

}elsif($fast =~ /N/i and $bam){
	$transformtobed_bam = "-b $bam";
}elsif($fast =~ /Y/i){
	$transformtobed_bam = "" ; # run transform_to_bed/pl withot bam file provided
}



##### firstly extracting reads aligned at TE#####

$cmd = "perl $bindir/lean_fq.pl -1 $rs1_ori -2 $rs2_ori -p $tmp_dir/rds_te -i $tmp_dir/$te_base -c $cpu_bwa "; 
process_cmd($cmd);
my $rds = "$tmp_dir/rds_te.fq1 $tmp_dir/rds_te.fq2";


#######  align reads associate with TE to the merged reference sequence  #######

$cmd = "bwa mem -T 20 -v 1 -t $cpu_bwa $tmp_dir/$proj.ref_and_te.fa $rds  |tee  $tmp_dir/$proj.ref_and_te.sam | samtools view -@ $cpu_view -bS - | samtools sort -@ $cpu_sort  - $tmp_dir/$proj.ref_and_te.sorted";


if (-e "$tmp_dir/$proj.ref_and_te.sam"){
	print STDERR "Seems alignment file exists. Skipped\n";
}else{
	process_cmd($cmd);
}

$cmd = "samtools index $tmp_dir/$proj.ref_and_te.sorted.bam";
process_cmd($cmd);

###### check  te seq id in te seq file####

open TE,$te_seq or die $!;
my @tes;
while(<TE>){
	chomp;
	if($_ =~ /^>(\S+)/){
		push @tes,$1;
	}
}

foreach my $te (@tes){

	########## extract informative reads from sam file  ###############

	$cmd = "perl $bindir/extract_informative.pl -s $tmp_dir/$proj.ref_and_te.sam  -n $te -p $tmp_dir/$proj  "; 
	process_cmd($cmd);
	$cmd = "perl $bindir/modify_informative.pl $tmp_dir/$proj.$te.informative.sam > $tmp_dir/$proj.$te.informative.full.sam";
	process_cmd($cmd);
	$cmd = "samtools view -bS $tmp_dir/$proj.$te.informative.sam | samtools sort - $tmp_dir/$proj.$te.informative.sorted";
	process_cmd($cmd);
	$cmd = "samtools index $tmp_dir/$proj.$te.informative.sorted.bam";
	process_cmd($cmd);


	#######   raln the informative reads back to the TE sequence  ###########
	$cmd = "perl $bindir/te_realin_bwa.pl -n $te -s $tmp_dir/$proj.$te.informative.full.sam -i $tmp_dir/$te_base -p $tmp_dir/$proj";
	process_cmd($cmd);


	######  identify the reads support insertion  #####
	$cmd = "perl $bindir/identity_inser_home_tnt_aa_debt.pl -s $tmp_dir/$proj.$te.informative.full.sam -g $tmp_dir/$proj.ref_and_te.fa -l $lib_len -n $te -r $tmp_dir/$proj.$te.alnte.sam -p $tmp_dir/$proj -a $lost";
	process_cmd($cmd);

	$cmd = "samtools view -bS $tmp_dir/${proj}.$te.supported.reads.sam | samtools sort - $tmp_dir/${proj}.$te.supported.reads.sorted ";
	process_cmd($cmd);
	$cmd = " samtools index   $tmp_dir/${proj}.$te.supported.reads.sorted.bam " ;
	process_cmd($cmd);
	###### 

	###### sort the support reads and generate bed files  ######
	$cmd = " sort -k 3,3 -k 4,4n $tmp_dir/${proj}.$te.ins.loc.lst >$tmp_dir/${proj}.$te.ins.loc.sorted.lst";
	process_cmd($cmd);
	$cmd = "perl  $bindir/transform_to_bed.pl $transformtobed_bam -n $te -p $tmp_dir/$proj.$te  -i $tmp_dir/$proj.$te.ins.loc.sorted.lst -w $window ";
	process_cmd($cmd);

	$cmd = "perl $bindir/filter_insertion.pl $para_filter -i $tmp_dir/$proj.$te.raw.bed -n $min_reads -q $map_q -r $ratio -d $depth_range >$tmp_dir/$proj.$te.filtered.bed";
	process_cmd($cmd);

	#####  Intergrate gene information in GFF and generate IGV snpshot batch file  ####
	if($gff){
		$cmd = "perl $bindir/annotate_bed.pl -b $tmp_dir/$proj.$te.filtered.bed  -a $gff -g $tmp_dir/$proj.ref_and_te.fa  -n $te -p $proj  -d $tmp_dir "; 
		process_cmd($cmd);
	}
}

######


	
	
#  subfuctions derived  from trinity package
sub process_cmd {
    my ($cmd) = @_;

    if($only_cmd =~ /Y/i){
		print STDERR "$cmd\n\n";
	}else{
		print STDERR &mytime."CMD: $cmd\n";
		my $start_time = time();
		my $ret = system($cmd);
		my $end_time = time();

		if ($ret) {
			die "Error, cmd: $cmd died with ret $ret";
		}
    
		print STDERR "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

		return;
	}
}





sub mytime() {
  my @mabbr = qw(January February March April May June July August September October November December);
  my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
  my $sec = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
  my $min = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
  my $hour = localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
  my $wday = $wabbr[localtime->wday];
  my $mday = localtime->mday;
  my $mon = $mabbr[localtime->mon];
  my $year = localtime->year() + 1900;
  return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}
