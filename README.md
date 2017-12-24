#ITIS:(Identify Transposon Insertion Sites)<br>
is a pipeline to identify novel TE insertion sites in genome

It require three input files:<br>  	
	(i) reference genome sequence,<br> 
	(ii)TE sequence, <br>
	(iii)paired-end (PE) short reads, with no restriction on length,generated from the re-sequenced genome that contains novel TE insertions.<br> 

By aligning read pairs to merged reference sequence, reference genome and TE sequence, ITIS will check each informative read pairs as long as it have more than 20bp overlap with TE sequence and determine if it supports the TE insertion around the location mapped by one of read pair.  In theory, by inspecting both cross read pairs and clipped reads at the same time, ITIS will have a higher sensitivity than other tools


NOTE: ITIS should be used to identy de novo insertions sites, It is unable to detect the lose event of one preexisting TEs.  
---
---


##Table of Contents
### <a href="#dep">Dependencies</a><br>
### <a href="#cmd">Command line options</a>
### <a href="#qck">Quick start with a sample dataset</a>
### <a href="#iss">Report an Issue</a>

---
---
### <a name="dep">Dependencies:

   The following programs need to be installed and the executable commands should be in $PATH of system:
	
	samtools (v 0.1.19)   #####******IMPORTANT******######
	bwa      (v 0.7.7-r441)
	bedtools (v 2.17.0)
	Bio::Perl
	blast+
	R
	
	Other usefull tool:
	IGV
	
-------------

### <a name="cmd">Comamnd Line Options

--------------------
USAGE:
#### perl itis.pl	
        /psc/home/jiangchuan/Dropbox/Code/Code_TE_inser/itis.pl
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
                        -d <Num,Num> the depth range to filter raw insertion site, [default 2,200]

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


                eg: perl  /psc/home/jiangchuan/Dropbox/Code/Code_TE_inser/itis.pl -g genome.fa -t tnt1.fa -l 300 -N test_run -1 reads.fq1 -2 reads.fq2 -f medicago.gff3 

-------------


### <a name="#qck">Quick start with a sample dataset
this test dataset derived from the genome resequencing project of Japonica A123(SRR631734), which have transposon mping be activted.

All PE reads mapped at chr1:1-2000000 were extracted and saved in file sample.fq1 and sample.fq2. 

First of all, untar the sample dataset:       
    
    cd test_dir
    tar xvzf sample_data.tar.gz     

####Input Files  
	PE reads:
		sample.fq1
		sample.fq1
	reference genome:
		rice_chr1_200k.fa
	Transposon sequence:
		mping.fa

####command to detect mping insertions in reference genome
	perl path/to/itis.pl -g rice_chr1_200k.fa -t mping.fa  -l 500 -N test -1 sample.fq1  -2 sample.fq2 -e Y    

	\#-e Y : to tell itis.pl that there are mping homologous sequence in reference genome

####Output Files
	itis will produce a lot of files in a directory named test.[aStringOfNumbers]
	
	The important files included:

		test.mping.filtered.bed  
			This is a list of reliable insertion sites.
			TAGS in column 4: 
				SR=(NO. of supporting library fragment),(No. supporting reads),(No. clipped reads at TE start),(No. clipped reads at TE end),(No. cross reads at TE start),(No. cross reads at TE end)
				MQ=average mapping quality of all supporting reads
				NM=name of TE
				GT=(No. of supporting reads),(No. of background reads):(Heter|Homo)
				PV=the p-value ofbinomial test for zygosity based on the GT values
				DP=average depth of 100bp region flanking the TE insertion site
				TS=the joint position at the begin of TE
				TE=the joint pisition at the end of TE
		                 result contents with test_data:
                     Chr1:0-2000000  174497  174500  SR=8,9,6,2,0,1;MQ=57;NM=mping;GT=8,0:Heter;PV=0.03125;DP=16;TS=1;TE=430;NB=N    .       +
                     Chr1:0-2000000  214352  214355  SR=15,17,9,5,1,2;MQ=57;NM=mping;GT=15,0:Homo;PV=0.00390625;DP=17;TS=1;TE=430;NB=N       .       +
                     Chr1:0-2000000  316534  316537  SR=9,14,3,11,0,0;MQ=57;NM=mping;GT=9,0:Heter;PV=0.03125;DP=20;TS=1;TE=430;NB=N  .       -
                     Chr1:0-2000000  639972  639975  SR=8,8,4,3,0,1;MQ=51;NM=mping;GT=8,14:Heter;PV=0.9903946;DP=28;TS=1;TE=430;NB=N .       +
                     Chr1:0-2000000  1193504 1193507 SR=9,12,5,7,0,0;MQ=56;NM=mping;GT=9,0:Heter;PV=0.03125;DP=18;TS=1;TE=430;NB=N   .       +
                     Chr1:0-2000000  1374936 1374939 SR=8,10,7,2,0,1;MQ=59;NM=mping;GT=8,10:Heter;PV=0.9407654;DP=28;TS=1;TE=430;NB=N        .       +
		test.mping.raw.bed  
			This is all candidate insertion sites, some of which may be false 
		test.mping.support.reads.sam and test.mping.support.reads.sorted.bam
			This is alignment file of all supportive reads
		test.all_reads_aln_ref_and_te.sort.bam
			This is alignment file of all reads
		commands_rcd
			A record of all the command used by itis.pl to identify TE insertions.  
		test.ref_and_te.fa    
			the reference sequence, containing genome and mping sequence,  used by bwa to align reads.   
	    
    *bam, *bed and reference sequence can be visuallized in IGV
	If you want to filter the raw insertion list by personalized criteria, you can rerun the script filter_insertion.pl, just as shown in command_rcd


	
-------------

### <a name="#iss">Report an Issue
If you have any questions or suggestion, please feel free to contact me :chuan-j@foxmail.com  or chjiang at sibs.ac.cn
-----------




