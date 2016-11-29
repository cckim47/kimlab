#!/usr/bin/perl -w

use strict;

my $indexKallisto = '/opt/index/kallisto/hs/genes.idx';
my $kallistodir = 'KaKa4Su';
my $kallistosam = 'Ka';
my $cpu = 4;

opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$kallistodir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $kallistodir;
    mkdir $kallistosam;
    system("kallisto quant -i $indexKallisto -t $cpu -b 100 -o $kallistodir --pseudobam --single -l 360 -s 200 $read1 > $kallistosam/aligned.sam");

    chdir('..');
}
