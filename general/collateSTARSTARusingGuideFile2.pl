#!/usr/bin/perl -w

use strict;

die "Usage: collateSTARSTARusingGuideFile.pl <guide_file>\n" unless @ARGV == 1;

open(GUIDE,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
my %files;
my @sampleOrder;
while(<GUIDE>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(' ');
    $files{$line[1]} = $line[0];
    push @sampleOrder,$line[1];
}
close GUIDE;

my %samplesGene;
my %geneCounts;
foreach my $sample (@sampleOrder) {
    $samplesGene{$sample}++;
    open(IN,$files{$sample}) or die "Couldn't open $files{$sample}: $!\n";
    chomp($_=<IN>); #header
    chomp($_=<IN>); #header
    chomp($_=<IN>); #header
    chomp($_=<IN>); #header
    while(<IN>) {
	chomp;
	$_ =~ s/\r//;
	next if !$_;
	my @line = split(/\t/);
	$geneCounts{$line[0]}{$sample} = $line[1];
#	    print "$line[0]\t$line[1]\n";
    }
    close IN;
}

# genes, Counts to stdout

my $sampleCount = 0;
foreach my $sampleID (@sampleOrder) {
    print "\t$sampleID";
    $sampleCount++;
}
print "\n";

foreach my $id (sort keys %geneCounts) {
    print "$id";
    foreach my $sampleID (@sampleOrder) {
        print "\t$geneCounts{$id}{$sampleID}";
    }
    print "\n";
}
