#!/usr/bin/perl -w

use strict;

my $indexSalmonQuasi = '/opt/index/salmon/hs/genes.quasi.index';
my $salmonQuasidir = 'SqSn4Su';
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

    chdir('..');
}
