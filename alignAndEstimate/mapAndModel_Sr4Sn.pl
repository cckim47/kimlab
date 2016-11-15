#!/usr/bin/perl -w

use strict;

my $genomeSTAR = '/home/charlie/index/star/hs';
my $outdir = 'Sr4Sn';
my $cpu = 12;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$outdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $outdir;
    system("STAR --runThreadN $cpu --outFileNamePrefix $outdir/ --genomeDir $genomeSTAR --readFilesIn $read1 --quantMode TranscriptomeSAM");

    chdir('..');
}

system("notify2.pl $outdir\_MappingComplete\n");
