#!/usr/bin/perl -w

use strict;

my $genesBt2 = '/opt/index/bowtie2_rseq/hs/genes_longName_bwt2.fa';
my $bwdir = 'Bw4Rs';


opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$bwdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $bwdir;
    system("bowtie2 -k 60 -i S,1,1.25 --score-min C,-14,0 --no-hd -x $genesBt2 -U $read1 -S $bwdir/mapped.reads.sam");
 
    chdir('..');
}
