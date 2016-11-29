#!/usr/bin/perl -w

use strict;

my $genomeRsem = '/opt/index/rsem/hs/genome';
my $outdir = 'SrRm';
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
    system("rsem-calculate-expression --star --star-path=/opt/bin/STAR -p $cpu $read1 $genomeRsem $subdir");
    system("mv *.results $outdir");
    system("mv *.transcript.* $outdir");
    system("mv *.stat $outdir");

    chdir('..');
}
