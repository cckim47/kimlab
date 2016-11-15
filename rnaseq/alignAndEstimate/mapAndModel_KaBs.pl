#!/usr/bin/perl -w

use strict;

my $kabsdir = 'KaBs';

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$kabsdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $kabsdir;
    system("Rscript ../getBitSeqExpression.R Ka");

    chdir('..');
}

system("notify2.pl KaBsMappingComplete\n");
