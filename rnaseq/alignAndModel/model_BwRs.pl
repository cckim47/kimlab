#!/usr/bin/perl -w

use strict;

my $genomeRseq = '/opt/index/hs/genes_longName.fa';
my $inDir = 'Bw4Rs';
my $outDir = 'BwRs';

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$outDir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $outDir;
    system("rseq expression_analysis $genomeRseq $indir/mapped.reads.sam");
    system("mv $indir/mapped.reads.sam.* $outDir/");
 
    chdir('..');
}
