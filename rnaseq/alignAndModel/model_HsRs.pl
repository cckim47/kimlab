#!/usr/bin/perl -w

use strict;

my $genomeRseq = '/opt/index/hs/genes_longName.fa';
my $indir = 'Hs4Rs';
my $infile = 'aligned.sam';
my $outdir = 'HsRs';

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
    system("rseq expression_analysis $genomeRseq $indir/$infile");
    system("mv $indir/$infile.* $outdir/");
 
    chdir('..');
}
