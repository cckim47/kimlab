#!/usr/bin/perl -w

use strict;

my $genes = '/home/charlie/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $thdir = 'Th';
my $Cudir = 'ThCu';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$Cudir/");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $Cudir;
    system("cufflinks -p $cpu -G $genesgtf -o $Cudir $thdir/accepted_hits.bam");

    chdir('..');
}
