#!/usr/bin/perl -w

use strict;

my $genome = '../../index/mm_grcm38/genome';
my $genomeRsem = '../../index/rsem/mm_grcm38/genome';
my $genes = '../../index/mm_grcm38/genes';
my $genesgtf = $genes . '.gtf';
my $genomeSTAR = '../../index/star/mm_grcm38';
my $thdir = '1_th_greedy';
my $cuffgdir = '2_cuffg';
my $cuffGdir = '3_cuffG';
my $htdir = '4_htseq';
my $rsemdir = '5_rsem';
my $stardir = '6_star';
my $cpu = 12;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$stardir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);
#    my $read2 = $subdir . '_R2' . ".fastq";

#    mkdir("$thdir");
#    system("tophat -p $cpu --transcriptome-index=$genes -o $thdir $genome $read1");

#    mkdir $cuffgdir;
#    system("cufflinks -p $cpu -g $genesgtf -o $cuffgdir $thdir/accepted_hits.bam");

#    mkdir $cuffGdir;
#    system("cufflinks -p $cpu -G $genesgtf -o $cuffGdir $thdir/accepted_hits.bam");

#    mkdir $htdir;
#    system("samtools view -h $thdir/accepted_hits.bam -o $htdir/accepted_hits.sam");
#    system("htseq-count --stranded=no $htdir/accepted_hits.sam $genesgtf > $htdir/accepted_hits.counts");

#    mkdir $rsemdir;
#    system("rsem-calculate-expression --bowtie2 -p $cpu $read1 $genomeRsem $subdir");
#    system("mv *.results $rsemdir");
#    system("mv *.transcript.* $rsemdir");
#    system("mv *.stat $rsemdir");

    mkdir $stardir;
    system("STAR --runThreadN $cpu --outFileNamePrefix $stardir/ --genomeDir $genomeSTAR --readFilesIn $read1 --quantMode GeneCounts");

    chdir('..');
}

system("notify2.pl STARmappingComplete\n");
