#!/usr/bin/perl -w

use strict;

my $indexSalmonFMD = '/opt/index/salmon/hs/genes.fmd.index';
my $salmonFMDdir = 'SfSn4Su';
my $cpu = 12;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$salmonFMDdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $salmonFMDdir;
    system("salmon quant -i $indexSalmonFMD --numBootstraps 100 -l U -p $cpu -r $read1 -o $salmonFMDdir");


    chdir('..');
}
