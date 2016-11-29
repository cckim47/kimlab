#!/usr/bin/perl -w

use strict;

my $genesgtf = '/opt/index/hs/genes.gtf';
my $indir = 'Th';
my $infile = 'accepted_hits.bam';
my $outdir = 'ThSt';
my $outfile = 'expression.gtf';
my $tableOut = 'expression.tab';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$outdir/");
    next unless (-e "$subdir/$subdir\.fastq");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $outdir;
    system("stringtie $indir/$infile -p $cpu -o $outdir/$outfile -G $genesgtf -A $outdir/$tableOut -B -e");

    chdir('..');
}
