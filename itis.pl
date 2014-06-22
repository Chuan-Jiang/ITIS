#!/usr/bin/perl
use warnings; use strict;
use FindBin;
use Getopt::Std;
use File::Basename;
use Time::localtime;

#########################  parameters #####################
my $usage = "USAGE:
	$0
	-g the genome sequence in fasta format  Required
	-t the TE sequence in fasta format		Required
	-T the temperate directory or use the defualt one	
	-n the id of your TE seq 	Required
	-N the name of you project Required
	-r the paired reads files seperate by space and quoted with ''    Required
	-h print STDERR this help message    
	-c cpu number for BWA   defualt 8i
	-l the average length of fragments in library	
	-w window size for cluster you support reads: default:50bp


	BWA should in you PATH
	";

die "$usage\n" if (@ARGV == 0);
my %opt;
getopts("g:t:T:n:N:r:c:l:h",\%opt);

die "$usage\n" if ($opt{h});

########===============================

my $genome = $opt{g};
my $te_seq = $opt{t};
my $te_base = basename $te_seq;
my $te     = $opt{n};
my $rs     = $opt{r};
my $cpu    = $opt{c}?$opt{c}:8;
my $proj   = $opt{N};
my $bindir  = "$FindBin::Bin";
my $lib_len = $opt{l};
my $tmp_dir = $opt{T}? $opt{T} : "tmp_".time();
my $window = $opt{w}?$opt{w}:50;
my $cmd;

##########################################################




####################### prepare template #########
mkdir $tmp_dir;					# temperate folder contain intermidiate files
$cmd = "cat $genome $te_seq >$tmp_dir/$te.genome.fa";
if (-e "$tmp_dir/$te.genome.fa"){
	print STDERR "Seems you already merged the sequences. Skipped\n";
}else{
	process_cmd($cmd);				# cat sequence together
}


$cmd = "cp $te_seq $tmp_dir/";
if (-e "$tmp_dir/$te_base"){
	print STDERR "Seems TE sequence already in the temperate directory, Skiped\n";
}else{
	process_cmd($cmd);				# copy te sequence to tmp/
}

###################### Index you sequence file #######
$cmd = "bwa index $tmp_dir/${te}.genome.fa";
if ( -e "$tmp_dir/$te.genome.fa.bwt"){
	print STDERR "Seems the Indexes for merged sequence exists. Skipped\n";
}else{
	process_cmd($cmd);				# index merged sequence
}

$cmd = "bwa index $tmp_dir/$te_base";
if ( -e "$tmp_dir/$te_base.bwt"){
	print STDERR "Seems the Indexes for TE sequence exists. Skipped\n";
}else{
	process_cmd($cmd);				# index te sequence
}


#######  align reads to the merged sequence and extract informative reads #######
$cmd = "bwa mem -v 1 -t $cpu $tmp_dir/${te}.genome.fa $rs > $tmp_dir/${te}.genome.sam";
process_cmd($cmd);
$cmd = "perl $bindir/extract_informative.pl -s $tmp_dir/${te}.genome.sam  -n $te -o $tmp_dir/$proj  "; 
process_cmd($cmd);
$cmd = "perl $bindir/modify_informative.pl $tmp_dir/$proj.informative.sam > $tmp_dir/$proj.informative.full.sam";
process_cmd($cmd);

#######   raln the informative reads back to the TE sequence  ###########
$cmd = "perl $bindir/te_realin_bwa.pl -t $te -s $tmp_dir/$proj.informative.sam -i $tmp_dir/$te_base -p $tmp_dir/$proj";
process_cmd($cmd);


######  identify the reads support insertion  #####
$cmd = "perl $bindir/identity_inser_home_tnt_aa_debt.pl -s $tmp_dir/$proj.informative.full.sam -g $tmp_dir/$te.genome.fa -i $lib_len -t $te -r $tmp_dir/$proj.alnte.sam -p $tmp_dir/$proj";
process_cmd($cmd);


###### sort the support reads and generate bed files  ######
$cmd = " sort -k 3,3 -k 4,4n $tmp_dir/${proj}_ins_loc.lst >${proj}_ins_loc_sorted.lst";
process_cmd($cmd);
$cmd = "perl  $bindir/transform_to_bed.pl  ${proj}_ins_loc_sorted.lst $window";
process_cmd($cmd);
$cmd = "samtools view -bS $tmp_dir/${proj}_supported_reads.sam >$tmp_dir/${proj}_supported_reads.bam";
process_cmd($cmd);
$cmd = "samtools sort $tmp_dir/${proj}_supported_reads.bam  ${proj}_supported_reads_sorted" ;
process_cmd($cmd);

######























	
	
sub process_cmd {
    my ($cmd) = @_;

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
