#!/usr/bin/perl -w

use strict;

my $indexSailfish = '/opt/index/sailfish/hs/';
my $sailfishdir = 'SlSl4Su';
my $cpu = 8;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$sailfishdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $sailfishdir;
   system("sailfish quant -i $indexSailfish --numBootstraps 100 -l U -p $cpu -r $read1 -o $sailfishdir");


    chdir('..');
}
