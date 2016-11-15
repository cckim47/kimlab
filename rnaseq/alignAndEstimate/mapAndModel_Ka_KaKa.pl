#!/usr/bin/perl -w

use strict;

my $genome = '/home/charlie/index/hs/genome';
my $genomeRsem = '/home/charlie/index/rsem/hs/genome';
my $genes = '/home/charlie/index/hs/genes';
my $genesgtf = $genes . '.gtf';
my $indexKallisto = '/opt/index/kallisto/hs/genes.idx';
my $kallistodir = 'KaKa';
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
#    system("kallisto quant -i $indexKallisto -o $kallistodir -t $cpu --single -l 360 -s 200 $read1");
    system("kallisto quant -i $indexKallisto -o $kallistodir --pseudobam --single -l 360 -s 200 $read1 > $kallistosam/aligned.sam");

    chdir('..');
}

system("notify2.pl kallistoMappingComplete\n");
