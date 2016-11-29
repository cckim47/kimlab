#!/usr/bin/perl -w

use strict;

## Need to convert mapping to transcripts (ENST##) to a format accepted by RSeq
## To do this, first gather gene/transcript mapping from  a conversion table listing transcripts and genes (conversion_tg.txt)
## Then, replace all instances of ENST## in the sam file to ENSG##$$ENST##

my $indir = 'Sr4Sn';
my $infile = 'Aligned.toTranscriptome.out.bam';
my $outdir = 'Sr4Rs';
my $outfile = 'Aligned.toTranscriptome.out.sam';
my $conversion = '../conversion_tg.txt';

# Read in conversion information
my %tid2longid;
open(ANN,$conversion) or die;
chomp($_=<ANN>); # header
while(<ANN>) {
    chomp;
    my @line = split(/\t/);
    $tid2longid{$line[0]} = $line[1] . '$$' . $line[0];
}
close ANN;

## Replace transcript names in each sam file
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

    # first convert bam to sam
    system("samtools view -h -o $outdir/$outfile $indir/$infile");

    # change annotations
    my $outfile2 = $outfile;
    $outfile2 =~ s/\.sam$/\.4Rs\.sam/;
    
    open(IN,"$outdir/$outfile") or die;
    open(OUT,">$outdir/$outfile2") or die;
    while(<IN>) {
	chomp;
	my $line = $_;
	my @idsToReplace;
	while ($line =~ /(ENST\d{11})/g) {
	    push @idsToReplace, $1;
	}
	foreach (@idsToReplace) {
	    $line =~ s/$_/$tid2longid{$_}/;
	}
	print OUT "$line\n";
    }
    
    chdir('..');
}
