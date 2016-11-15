#!/usr/bin/perl -w

use strict;

my $genomeRsem = '/opt/index/rsem/hs/genome';
my $rsemstardir = 'SrRm';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$rsemstardir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);


    mkdir $rsemstardir;
    system("rsem-calculate-expression --star --star-path=/opt/bin/STAR -p $cpu $read1 $genomeRsem $subdir");
    system("mv *.results $rsemstardir");
    system("mv *.transcript.* $rsemstardir");
    system("mv *.stat $rsemstardir");

    chdir('..');
}

system("notify2.pl $rsemstardir\MappingComplete\n");
