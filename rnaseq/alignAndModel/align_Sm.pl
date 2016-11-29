#!/usr/bin/perl -w

use strict;

my $genomeSeqmap = '/opt/index/hs/genes_longName.fa';
my $seqmapDir = 'Sm'; 

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$seqmapDir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $seqmapDir;
    system("seqmap 2 $read1 $genomeSeqmap $seqmapDir/mapped.reads.out.sam /eland:3 /available_memory:12000");

    chdir('..');
}
