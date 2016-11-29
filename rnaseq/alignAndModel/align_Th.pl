#!/usr/bin/perl -w

use strict;

my $genome = '/home/charlie/index/hs/genome';
my $genes = '/home/charlie/index/hs/genes';
my $thdir = 'Th';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$thdir/");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir("$thdir");
    my $read1 = $subdir . ".fastq";
    system("tophat -p $cpu --transcriptome-index=$genes -o $thdir $genome $read1");

    chdir('..');
}
