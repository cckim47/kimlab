#!/usr/bin/perl -w

use strict;

my $genes = '/opt/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $stardir = 'Sr';
my $htSTARdir = 'SrHt';

my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$htSTARdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);
#    my $read1 = $subdir . ".fastq";

    mkdir $htSTARdir;
    system("htseq-count --stranded=no $stardir/Aligned.out.sam $genesgtf > $htSTARdir/Aligned.out.counts");

    chdir('..');
}

system("notify2.pl $htSTARdir\MappingComplete\n");
