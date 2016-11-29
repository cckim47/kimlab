#!/usr/bin/perl -w

use strict;

my $genomeXs = '/opt/index/hs/genes.fa';
my $indir = 'Hs';
my $infile = 'aligned.sam';
my $outdir = 'HsXs';

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
    system("express $genomeXs $indir/$infile -m 360 -s 200 -o $outdir");

    chdir('..');
}
