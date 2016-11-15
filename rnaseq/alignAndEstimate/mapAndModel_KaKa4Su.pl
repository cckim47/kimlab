#!/usr/bin/perl -w

use strict;

#my $genome = '/home/charlie/index/hs/genome';
#my $genomeRsem = '/home/charlie/index/rsem/hs/genome';
#my $genes = '/home/charlie/index/hs/genes';
#my $genesgtf = $genes . '.gtf';
my $indexKallisto = '/opt/index/kallisto/hs/genes.idx';
my $kallistodir = 'KaKa4Su';
my $cpu = 10;

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
    system("kallisto quant -i $indexKallisto -t $cpu -b 100 -o $kallistodir --single -l 360 -s 200 $read1");

    chdir('..');
}

system("notify2.pl KaKa4SuMappingComplete\n");
