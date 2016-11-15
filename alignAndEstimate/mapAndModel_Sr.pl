#!/usr/bin/perl -w

use strict;

my $genome = '/home/charlie/index/hs/genome';
my $genomeRsem = '/home/charlie/index/rsem/hs/genome';
my $genes = '/home/charlie/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $genomeSTAR = '/home/charlie/index/star/hs';
my $stardir = 'Sr';
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

    mkdir $stardir;
    system("STAR --runThreadN $cpu --outFileNamePrefix $stardir/ --genomeDir $genomeSTAR --readFilesIn $read1 --quantMode GeneCounts");

    chdir('..');
}

system("notify2.pl starMappingComplete\n");
