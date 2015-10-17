#!/usr/bin/env perl 
use warnings; use strict;
use FindBin;
use Getopt::Std;
use File::Basename;
use Time::localtime;

#########################  parameters #####################
my $usage = "USAGE:
	$0
	REQUIRED -g the genome sequence file in fasta format  
		OR -G prefix of bwa-indexed reference file ( genome + transposon) 
	REQUIRED -t the TE sequence file in fasta format
		or -T prefix of bwa-indexed transposon sequence file
	REQUIRED -l the average length of fragments in library	
	REQUIRED -N the name of project
	REQUIRED -1 the paired read file 1
	REQUIRED -2 the paired read file 2	
		
		-f <gff file> if provided, ITIS will check if TE inserted in gentic or intergeneic region
		-F <Y|N> run scripts in 'FAST' mode; It won't align all reads to reference genome,caculate the average bg depth,
		     and estimate if insertion is homo or heter,[default N]
			
			##  parameters used with  '-F N'  :
			-B <bam file> use previous sorted and indexed bam file of reads aligned to reference genome
			-D <Num,Num> the depth range to filter raw insertion site, [default 2,200]

		-q <Num>  the minimum average mapping quality of all supporting reads, [default 1]

		-e <Y|N> If reference genome contains this TE or it's homolog. using blast to hard mask these sequence is required, [default N]
		
		-a <Num> the number of bases allowed to be lost when transposing, [defualt 10]

		-b <tags> minimum required number of flanking reads , in the format of /Tag=Value/Tag=Value/Tag=Value/ , the avaliable tags:
			t: total supporting reads at detected insertion  /t=3/
			CS:clipped reads cover TE start site /CS=0/
			CE:clipped reads cover TE end site  /CE=0/
			cs:cross reads cover TE start  /cs=0/
			ce:cross reads cover TE end    /cs=0/
			TS:total reads cover TE start  /TS=1/
			TE:total reads cover TE end    /TE=1/
				[default /t=3/TS=1/TE=1/]
		-c <Num,Num,Num> cpu number for 'BWA mem', 'samtools view'  and 'samtools sort', [defualt 8,2,2]
		
		-w <Num> window size used to cluster supportting reads, [default library_length/2]
		
		-D <Directory> use this specifed temperate directory, [default[project].[aStringOfNumbers]]
		
		-m <Y|N> Only print out all commands to STDERR, [default N]
	         
		-h print this help message    

	
		eg: perl  $0 -g genome.fa -t tnt1.fa -l 300 -N test_run -1 reads.fq1 -2 reads.fq2 -f medicago.gff3 

	
	";

die "$usage\n" if (@ARGV == 0);
my %opt;
getopts("g:t:l:N:1:2:f:b:R:B:D:G:T:c:q:e:a:F:D:w:m:h",\%opt);

die "$usage\n" if ($opt{h});

########===============================

my $genome = $opt{g};
my $te_seq = $opt{t};

my $index_te = $opt{T};
my $index_ref = $opt{G};

my $lib_len = $opt{l};
my $proj   = $opt{N};
my $rs1_ori  = $opt{1};
my $rs2_ori  = $opt{2};
my $gff = $opt{f};

my $min_reads  = $opt{b}?$opt{b}:"/t=3/TS=1/TE=1/";
my $bam = $opt{B}?$opt{B}:0;
my $depth_range= $opt{D}?$opt{D}:"2,200";
my $cpu    = $opt{c}?$opt{c}:"8,2,2";
	my($cpu_bwa,$cpu_view,$cpu_sort) = split /,/,$cpu;
my $exists = $opt{e}?$opt{e}:"N";
my $fast   = $opt{F}?$opt{F}:"N";
my $tmp_dir = $opt{D}? $opt{D} : "$proj.".time();
my $window = $opt{w}?$opt{w}:$lib_len/2;
my $only_cmd = $opt{m}?$opt{m}:"N";
my $cmd;
my $bindir  = "$FindBin::Bin";
my $map_q = $opt{q}?$opt{q}:1;
my $lost = $opt{a}?$opt{a}:10;
my $bwa = "bwa";
##########################################################
  


if(-e $tmp_dir){
	print "using dir: $tmp_dir\n";
}else{
	$cmd = "mkdir $tmp_dir";
	system($cmd) == 0 or die $!;
}

open CMD, ">$tmp_dir/commands_rcd" or die $!;

my $para_filter = "";

if($index_ref){
	$genome = $index_ref;
}else{
	####################### prepare reference  #########

	if($exists =~ /N/i){
		$cmd = "cat $genome $te_seq >$tmp_dir/$proj.ref_and_te.fa";
	}else{
		$cmd = "perl -I $bindir $bindir/mask_te_homo_in_genome.pl -g $genome -t $te_seq -p $tmp_dir/te_homo_in_ref.lst -o $tmp_dir/$proj.ref_and_te.fa";
		$para_filter = "-l $tmp_dir/te_homo_in_ref.lst";
	}

	process_cmd($cmd);				# cat sequence together
	
	
	###################### Index sequence file #######

	$cmd = "$bwa index $tmp_dir/$proj.ref_and_te.fa";
	if ( -e "$tmp_dir/$proj.ref_and_te.fa.bwt"){
		print STDERR "Seems like the Indexes for merged sequence exists. Skipped\n";
	}else{
		process_cmd($cmd);				# index merged sequence
	}

	$index_ref = "$tmp_dir/$proj.ref_and_te.fa";
}
if($index_te){
	1;
}else{
	$cmd = "cp $te_seq $tmp_dir/";
	my $te_base = basename $te_seq;
	process_cmd($cmd);              # copy te sequence to tmp/
	$cmd = "$bwa index $tmp_dir/$te_base";
	process_cmd($cmd);              # index te sequence
	$index_te = "$tmp_dir/$te_base";
}


##### align original reads to reference genome ######
my $transformtobed_bam ;
if($fast =~ /N/i and $bam == 0){
	$cmd = "$bwa mem -T 20 -t $cpu_bwa $index_ref $rs1_ori $rs2_ori 2>/dev/null | samtools view -@ $cpu_view -buS - | samtools sort -@ $cpu_sort - $tmp_dir/$proj.all_reads_aln_ref_and_te.sort";
	process_cmd($cmd);
	$transformtobed_bam = "-b $tmp_dir/$proj.all_reads_aln_ref_and_te.sort.bam";
	
	$cmd = "samtools index $tmp_dir/$proj.all_reads_aln_ref_and_te.sort.bam";
	process_cmd($cmd);

}elsif($fast =~ /N/i and $bam){
	$transformtobed_bam = "-b $bam";
}elsif($fast =~ /Y/i){
	$transformtobed_bam = "" ; # run transform_to_bed/pl withot bam file provided
}



##### firstly extracting reads aligned at TE#####

$cmd = "perl -I $bindir $bindir/lean_fq.pl -1 $rs1_ori -2 $rs2_ori -p $tmp_dir/rds_te -i $index_te -c $cpu_bwa "; 
process_cmd($cmd);
my $rds = "$tmp_dir/rds_te.fq1 $tmp_dir/rds_te.fq2";


#######  align reads associate with TE to the merged reference sequence  #######

$cmd = "$bwa mem -T 20 -v 1 -t $cpu_bwa $index_ref $rds 2>/dev/null  |tee  $tmp_dir/$proj.ref_and_te.sam | samtools view -@ $cpu_view -buS - | samtools sort -@ $cpu_sort  - $tmp_dir/$proj.ref_and_te.sorted";


if (-e "$tmp_dir/$proj.ref_and_te.sam"){
	print STDERR "Seems alignment file exists. Skipped\n";
}else{
	process_cmd($cmd);
}

$cmd = "samtools index $tmp_dir/$proj.ref_and_te.sorted.bam";
process_cmd($cmd);

###### check  te seq id in te seq file####

open TE,$index_te or die $!;

my @tes;
while(<TE>){
	chomp;
	if($_ =~ /^>(\S+)/){
		push @tes,$1;
	}
}

foreach my $te (@tes){

	########## extract informative reads from sam file  ###############

	$cmd = "perl -I $bindir $bindir/extract_informative.pl -g $genome -s $tmp_dir/$proj.ref_and_te.sam  -n $te -p $tmp_dir/$proj  "; 
	process_cmd($cmd);
	$cmd = "perl -I $bindir $bindir/modify_informative.pl $tmp_dir/$proj.$te.informative.sam > $tmp_dir/$proj.$te.informative.full.sam";
	process_cmd($cmd);
	$cmd = "samtools view -buS $tmp_dir/$proj.$te.informative.sam | samtools sort - $tmp_dir/$proj.$te.informative.sorted";
	process_cmd($cmd);
	$cmd = "samtools index $tmp_dir/$proj.$te.informative.sorted.bam";
	process_cmd($cmd);


	#######   raln the informative reads back to the TE sequence  ###########
	$cmd = "perl -I $bindir $bindir/te_realin_bwa.pl -n $te -s $tmp_dir/$proj.$te.informative.full.sam -i $index_te -p $tmp_dir/$proj";
	process_cmd($cmd);


	######  identify the reads support insertion  #####
	$cmd = "perl -I $bindir $bindir/identity_inser_sites.pl -s $tmp_dir/$proj.$te.informative.full.sam -g $index_ref -l $lib_len -n $te -r $tmp_dir/$proj.$te.alnte.sam -p $tmp_dir/$proj -a $lost";
	process_cmd($cmd);

	$cmd = "samtools view -buS $tmp_dir/${proj}.$te.support.reads.sam | samtools sort - $tmp_dir/${proj}.$te.support.reads.sorted ";
	process_cmd($cmd);
	$cmd = "samtools index   $tmp_dir/${proj}.$te.support.reads.sorted.bam " ;
	process_cmd($cmd);
	###### 

	###### sort the support reads and generate bed files  ######
	$cmd = "sort -k 3,3 -k 4,4n $tmp_dir/${proj}.$te.ins.loc.lst >$tmp_dir/${proj}.$te.ins.loc.sorted.lst";
	process_cmd($cmd);
	$cmd = "perl -I $bindir  $bindir/transform_to_bed.pl $transformtobed_bam -e $index_te -n $te -p $tmp_dir/$proj.$te  -i $tmp_dir/$proj.$te.ins.loc.sorted.lst -l $lib_len  -w $window ";
	process_cmd($cmd);

	$cmd = "perl -I $bindir $bindir/filter_insertion.pl $para_filter -i $tmp_dir/$proj.$te.raw.bed -n $min_reads -q $map_q  -d $depth_range >$tmp_dir/$proj.$te.filtered.bed";
	process_cmd($cmd);

	#####  Intergrate gene information in GFF and generate IGV snpshot batch file  ####
	if($gff){
		$cmd = "perl -I $bindir $bindir/annotate_bed.pl -b $tmp_dir/$proj.$te.filtered.bed  -a $gff -g $index_ref  -n $te -p $proj  -d $tmp_dir "; 
		process_cmd($cmd);
	}
}

######


	
	
#  subfuctions derived  from trinity package
sub process_cmd {
    my ($cmd) = @_;

    if($only_cmd =~ /Y/i){
		print STDERR "$cmd\n\n";
		print CMD " $cmd\n\n";
	}else{
		print STDERR &mytime."CMD: $cmd\n";
		print CMD "$cmd\t#".&mytime."\n";
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
