#!/usr/bin/perl -w

use strict;

## Calls an Rscript to run tximport on Kallisto, Salmon, or Sailfish transcript level estimates
## Specify the input directory (ex. KaKa4Su) in the command line for the script call

my $indir = $ARGV[0];
die "Usage: model_2g.pl <inputDir>\n" unless @ARGV == 1;
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
    system("Rscript ../runTximport.R $indir $gene");

    chdir('..');
}
