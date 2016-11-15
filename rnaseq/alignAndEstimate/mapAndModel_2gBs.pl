#!/usr/bin/perl -w

use strict;

my $indir = $ARGV[0];
die "Usage: mapAndModel_2gBs.pl <inputDir>\n" unless @ARGV == 1;
my $outdir = $indir . '2g';
my $gene = '../../conversion_tg.txt';

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
    system("Rscript ../runTximport_Bs.R $indir $gene");

    chdir('..');
}

system("notify2.pl $outdir\_MappingComplete\n");
