#!/usr/bin/perl -w

use strict;

my $genomeSeqmap = '/opt/index/hs/genes_longName.fa';
my $seqmapDir = 'Sm'; 
my $rseqDir = 'SmRs';

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$rseqDir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $rseqDir;
    ## Change the seqmap eland output into a format compatible with rseq
    system("tr ' ' '_' < $seqmapDir/mapped.reads.out.sam > $seqmapDir/mapped.reads.out_noSpace.sam");

    system("rseq expression_analysis $genomeSeqmap $seqmapDir/mapped.reads.out_noSpace.sam");
    system("mv $seqmapDir/mapped.reads.out_noSpace.* $rseqDir/");

    chdir('..');
}
