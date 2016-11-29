#!/usr/bin/perl -w

use strict;

my $genomeRsem = '/opt/index/rsem/hs/genome';
my $rsemdir = 'BwRm';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$rsemdir/");

    print "Analyzing $subdir\n";
    chdir($subdir);

    my $read1 = $subdir . ".fastq";
    mkdir $rsemdir;
    system("rsem-calculate-expression --bowtie2 -p $cpu $read1 $genomeRsem $subdir");
    system("mv *.results $rsemdir");
    system("mv *.transcript.* $rsemdir");
    system("mv *.stat $rsemdir");

    chdir('..');
}
