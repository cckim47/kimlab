#!/usr/bin/perl -w

use strict;

my $genes = '/home/charlie/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $thdir = 'Th';
my $htdir = 'ThHt';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$htdir/");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $htdir;
    ## Convert bam to sam
    system("samtools view -h $thdir/accepted_hits.bam -o $htdir/accepted_hits.sam");
    system("htseq-count --stranded=no $htdir/accepted_hits.sam $genesgtf > $htdir/accepted_hits.counts");

    chdir('..');
}
