#!/usr/bin/perl -w

use strict;

my $genomeSTAR = '/home/charlie/index/star/hs';
my $stardir = 'SrSr';
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
