#!/usr/bin/perl -w

use strict;

my $genesgtf = $genes . '.gtf';
my $indir = 'SrSr';
my $outdir = 'SrHt';

my $cpu = 4;

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
    system("htseq-count --stranded=no $indir/Aligned.out.sam $genesgtf > $outdir/Aligned.out.counts");

    chdir('..');
}
