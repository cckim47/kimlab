#!/usr/bin/perl -w

use strict;

my $genesgtf = '/opt/index/hs/genes.gtf';
my $indir = 'Hs4St';
my $infile = 'aligned.sam';
my $outdir = 'HsSt';
my $outfile = 'expression.gtf';
my $tableOut = 'expression.tab';
my $cpu = 2;

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
    system("samtools sort -o $outdir/$infile\.sorted $indir/$infile");
    system("stringtie $outdir/$infile\.sorted -p $cpu -o $outdir/$outfile -G $genesgtf -A $outdir/$tableOut -B -e");

    chdir('..');
}

system("notify2.pl $outdir\_MappingComplete\n");
