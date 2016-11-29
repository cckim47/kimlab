#!/usr/bin/perl -w

use strict;

my $genomeHisat2 = '/opt/index/hisat2_genome/hs/genome';
my $splices = '/opt/index/hisat2_genome/hs/splicesites.txt';
my $outdir = 'Hs4St';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$outdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $outdir;
    system("hisat2 --dta -p $cpu -x $genomeHisat2 -U $read1 --known-splicesite-infile $splices -S $outdir/aligned.sam");
 
    chdir('..');
}
