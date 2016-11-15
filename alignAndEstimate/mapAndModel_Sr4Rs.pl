#!/usr/bin/perl -w

use strict;

my $indir = 'Sr4Sn';
my $infile = 'Aligned.toTranscriptome.out.bam';
my $outdir = 'Sr4Rs';
my $outfile = 'Aligned.toTranscriptome.out.sam';
my $conversion = '../conversion_tg.txt';


my %tid2longid;
open(ANN,$conversion) or die;
chomp($_=<ANN>); # header
while(<ANN>) {
    chomp;
    my @line = split(/\t/);
    $tid2longid{$line[0]} = $line[1] . '$$' . $line[0];
#    print "$line[0]\t$tid2longid{$line[0]}\n";
}
close ANN;


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

    # first convert to sam
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

system("notify2.pl $outdir\_MappingComplete\n");
