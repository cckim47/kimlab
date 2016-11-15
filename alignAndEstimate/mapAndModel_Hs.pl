#!/usr/bin/perl -w

use strict;

my $genome = '/home/charlie/index/hs/genome';
my $genomeRsem = '/home/charlie/index/rsem/hs/genome';
my $genes = '/home/charlie/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $genomeHisat2 = '/opt/index/hisat2/hs/genes';
my $hisat2dir = 'Hs';


opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$hisat2dir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $hisat2dir;
    system("hisat2 -x $genomeHisat2 -U $read1 -S $hisat2dir/aligned.sam");
 
    chdir('..');
}

system("notify2.pl hisat2MappingComplete\n");
