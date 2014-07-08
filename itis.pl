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
	REQUIRED -n the seq id of your TE seq, must be identical to the one of TE seq file
	REQUIRED -N the name of you project
	REQUIRED -1 the paired read file 1
	REQUIRED -2 the paired read file 2	
	REQUIRED -f gff file
			 
		-b the total number of required events[3] OR minimum required support events at TE start[1], events at TE end[1], cross reads at TE start
			defualt in the form of '3,1,1', seperate by comma 
		-R the minimum ratio of support fragments  and depth of 200 bp around the insertion site; default 0.2
		-B the sorted and indexed bam file of all clean reads align to reference genome; on condition of '-F F'
		-D <3,200>, the depth range to filter candidate insertion site. 
		-c FORMAT:\\d,\\d,\\d; for  cpu number for BWA and samtools view and samtools sort    defualt 8,2,2
		-e <T|F: default F> if TE sequence have homolog in genome. use blast to hard mask repeat sequence
		-F <T|F: default F> run scripts in 'FAST' mode; It won't align all reads to reference genome , you can do it by yourself
		-T use this specifed temperate directory or use the DEFAULT one :[project].[aStringOfNumbers]
		-w window size for cluster you support reads: DEFAULT :100bp
		-m <T|F: default F> Only print out all commands to STDERR
	         
		-h print STDERR this help message    

	BWA samtools should in you PATH
	";

die "$usage\n" if (@ARGV == 0);
my %opt;
getopts("g:t:l:n:N:1:2:f:b:R:B:D:c:e:F:T:w:m:h",\%opt);

die "$usage\n" if ($opt{h});

########===============================

my $genome = $opt{g};
my $te_seq = $opt{t};
	my $te_base = basename $te_seq;
my $lib_len = $opt{l};
my $te     = $opt{n};
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
my $exists = $opt{e}?$opt{e}:"F";
my $fast   = $opt{F}?$opt{F}:"F";
my $tmp_dir = $opt{T}? $opt{T} : "$proj.".time();
my $window = $opt{w}?$opt{w}:100;
my $only_cmd = $opt{m}?$opt{m}:"F";
my $cmd;
my $bindir  = "$FindBin::Bin";


##########################################################
  

$cmd = "mkdir $tmp_dir";
process_cmd($cmd);					# temperate folder contain intermidiate files


####################### prepare template #########

if($exists =~ /F/i){
	$cmd = "cat $genome $te_seq >$tmp_dir/$proj.$te.genome.fa";
}else{
	$cmd = "perl $bindir/mask_te_homo_in_genome.pl -g $genome -t $te_seq -o $tmp_dir/$proj.$te.genome.fa";
}

if (-e "$tmp_dir/$proj.$te.genome.fa"){
	print STDERR "Seems you already merged the sequences. Skipped\n";
}else{
	process_cmd($cmd);				# cat sequence together
}


$cmd = "cp $te_seq $tmp_dir/";
process_cmd($cmd);				# copy te sequence to tmp/


###################### Index you sequence file #######


$cmd = "bwa index $tmp_dir/$proj.$te.genome.fa";
if ( -e "$tmp_dir/$proj.$te.genome.fa.bwt"){
	print STDERR "Seems the Indexes for merged sequence exists. Skipped\n";
}else{
	process_cmd($cmd);				# index merged sequence
}


$cmd = "bwa index $tmp_dir/$te_base";
process_cmd($cmd);				# index te sequence


##### align original reads to reference genome ######
my $transformtobed_bam ;
if($fast =~ /F/i and $bam == 0){
	$cmd = "bwa mem -T 20 -t $cpu_bwa $tmp_dir/$proj.$te.genome.fa $rs1_ori $rs2_ori | samtools view -@ $cpu_view -bS - | samtools sort -@ $cpu_sort - $tmp_dir/$proj.aln.bg.ref.sort";
	process_cmd($cmd);
	$transformtobed_bam = "-D $depth_range -b $tmp_dir/$proj.aln.bg.ref.sort.bam";
	
	$cmd = "samtools index $tmp_dir/$proj.aln.bg.ref.sort.bam";
	process_cmd($cmd);

}elsif($fast =~ /F/i and $bam){
	$transformtobed_bam = "-D $depth_range -b $bam";
}elsif($fast =~ /T/i){
	$transformtobed_bam = "" ; # run transform_to_bed/pl withot bam file provided
}



##### firstly extracting reads aligned at TE#####

$cmd = "perl $bindir/lean_fq.pl -1 $rs1_ori -2 $rs2_ori -p $tmp_dir/rds_te -i $tmp_dir/$te_base -c $cpu_bwa "; 
process_cmd($cmd);
my $rds = "$tmp_dir/rds_te.fq1 $tmp_dir/rds_te.fq2";


#######  align reads associate with TE to the merged reference sequence  #######

$cmd = "bwa mem -T 20 -v 1 -t $cpu_bwa $tmp_dir/$proj.$te.genome.fa $rds  |tee  $tmp_dir/$proj.$te.genome.sam | samtools view -@ $cpu_view -bS - | samtools sort -@ $cpu_sort  - $tmp_dir/$proj.$te.genome.sorted";


if (-e "$tmp_dir/$proj.$te.genome.sam"){
	print STDERR "Seems alignment file exists. Skipped\n";
}else{
	process_cmd($cmd);
}

$cmd = "samtools index $tmp_dir/$proj.$te.genome.sorted.bam";
process_cmd($cmd);

########## extract informative reads from sam file  ###############

$cmd = "perl $bindir/extract_informative.pl -s $tmp_dir/$proj.$te.genome.sam  -n $te -p $tmp_dir/$proj  "; 
process_cmd($cmd);
$cmd = "perl $bindir/modify_informative.pl $tmp_dir/$proj.informative.sam > $tmp_dir/$proj.informative.full.sam";
process_cmd($cmd);
$cmd = "samtools view -bS $tmp_dir/$proj.informative.sam | samtools sort - $tmp_dir/$proj.informative.sorted";
process_cmd($cmd);
$cmd = "samtools index $tmp_dir/$proj.informative.sorted.bam";
process_cmd($cmd);


#######   raln the informative reads back to the TE sequence  ###########
$cmd = "perl $bindir/te_realin_bwa.pl -n $te -s $tmp_dir/$proj.informative.full.sam -i $tmp_dir/$te_base -p $tmp_dir/$proj";
process_cmd($cmd);


######  identify the reads support insertion  #####
$cmd = "perl $bindir/identity_inser_home_tnt_aa_debt.pl -s $tmp_dir/$proj.informative.full.sam -g $tmp_dir/$proj.$te.genome.fa -l $lib_len -n $te -r $tmp_dir/$proj.alnte.sam -p $tmp_dir/$proj";
process_cmd($cmd);

$cmd = "samtools view -bS $tmp_dir/${proj}.supported.reads.sam | samtools sort - $tmp_dir/${proj}.supported.reads.sorted ";
process_cmd($cmd);
$cmd = " samtools index   $tmp_dir/${proj}.supported.reads.sorted.bam " ;
process_cmd($cmd);
###### 

###### sort the support reads and generate bed files  ######
$cmd = " sort -k 3,3 -k 4,4n $tmp_dir/${proj}.ins.loc.lst >$tmp_dir/${proj}.ins.loc.sorted.lst";
process_cmd($cmd);
$cmd = "perl  $bindir/transform_to_bed.pl -p $tmp_dir/$proj $transformtobed_bam -t $min_reads -r $ratio -i $tmp_dir/$proj.ins.loc.sorted.lst -w $window ";
process_cmd($cmd);


#####  Intergrate gene information in GFF and generate IGV snpshot batch file  ####
if($gff){
	$cmd = "perl $bindir/annotate_bed.pl -b $tmp_dir/$proj.filtered.bed  -a $gff -g $tmp_dir/$proj.$te.genome.fa  -n $te -p $proj  -d $tmp_dir "; 
	process_cmd($cmd);
}

######


	
	
#  subfuctions derived  from trinity package
sub process_cmd {
    my ($cmd) = @_;

    if($only_cmd =~ /T/i){
		print STDERR "CMD: $cmd\n";
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
