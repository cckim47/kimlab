#!/usr/bin/perl -w

use strict;

my $genome = '/home/charlie/index/hs/genome';
my $genomeRsem = '/home/charlie/index/rsem/hs/genome';
my $genes = '/home/charlie/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $indexSalmonQuasi = '/opt/index/salmon/hs/genes.quasi.index';
my $indexSalmonFMD = '/opt/index/salmon/hs/genes.fmd.index';
my $salmonQuasidir = 'SqSn4Su';
my $salmonFMDdir = 'SfSn4Su';
my $cpu = 12;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$salmonQuasidir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $salmonQuasidir;
    system("salmon quant -i $indexSalmonQuasi --numBootstraps 100 -l U -p $cpu -r $read1 -o $salmonQuasidir");

    mkdir $salmonFMDdir;
    system("salmon quant -i $indexSalmonFMD --numBootstraps 100 -l U -p $cpu -r $read1 -o $salmonFMDdir");


    chdir('..');
}

system("notify2.pl salmonMappingComplete\n");
