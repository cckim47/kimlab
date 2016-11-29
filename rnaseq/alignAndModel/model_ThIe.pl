#!/usr/bin/perl -w

use strict;

my $genes = '/opt/index/isoem/hs/genes_reduced.gtf';
my $indir = 'Th';
my $infile = 'accepted_hits.bam';
my $outdir = 'ThIe';
my $outfile = 'accepted_hits.sam';


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
    ## Convert bam to sam
    system("samtools view -h -o $outdir/$outfile $indir/$infile");
    system("isoem2 -G $genes -m 360 -d 200 -O $outdir $outdir/$outfile");
 
    chdir('..');
}
