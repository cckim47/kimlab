#!/usr/bin/perl -w

use strict;

## Bitseq is called from an Rscript for each sample

my $hsbsdir = 'HsBs';

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$hsbsdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $hsbsdir;
    system("Rscript ../model_Bs.R Hs");

    chdir('..');
}
