#!/usr/bin/perl -w

use strict;

my $genome = '../../../index/hs/genome';
my $genomeRsem = '../../../index/rsem/hs/genome';
my $genes = '../../../index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $genomeSTAR = '../../../index/star/hs';
my $thdir = '1_th_greedy';
my $cuffgdir = '2_cuffg';
my $cuffGdir = '3_cuffG';
my $htdir = '4_htseq';
my $rsemdir = '5_rsem';
my $stardir = '6_star';
my $htSTARdir = '7_starhtseq';
my $rsemstardir = '8_starrsem';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$rsemstardir/");

    print "Analyzing $subdir\n";
    chdir($subdir);

    my $read1 = $subdir . ".fastq";

    mkdir $rsemstardir;
    system("~/bin/rsem-1.2.22/rsem-calculate-expression --star --star-path ~/bin -p $cpu $read1 $genomeRsem $subdir");
    system("mv *.results $rsemstardir");
    system("mv *.transcript.* $rsemstardir");
    system("mv *.stat $rsemstardir");

    chdir('..');
}

system("notify2.pl MappingComplete\n");
