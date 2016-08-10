#!/usr/bin/perl -w

use strict;

die "Usage: updateAnnotation.pl <input_pcl_or_cdt> <input_ann>\n" unless @ARGV == 2;

my %ann;
open(ANN,$ARGV[1]) or die "Couldn't open $ARGV[1]: $!\n";
while(<ANN>) {
    chomp;
    next if !$_;
    $_ =~ s/\r//;
    my @line = split(/\t/);
    $ann{$line[0]} = $line[1];
}
close ANN;

my $gidFlag = 0;
open(DATA,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
chomp($_=<DATA>);
my @headSplit = split(/\t/);
my $columns = $#headSplit;
seek DATA,0,0;

while(<DATA>) {
    chomp;
    next if !$_;
    $_ =~ s/\r//;

    if (/^GID/) {
        $gidFlag++;
        print "$_\n";
        next;
    }
    elsif (/^UNIQID/) {
        print "$_\n";
        next;
    }
    elsif (/^id/) {
        print "$_\n";
        next;
    }
    elsif (/^EWEIGHT/) {
        print "$_\n";
        next;
    }
    elsif (/^AID/) {
        print "$_\n";
        next;
    }
    else {
        my $origLine = $_;
        my @line = split(/\t/,$origLine);
        $line[1+$gidFlag] = $ann{$line[0+$gidFlag]};

	unless ($line[1+$gidFlag]) {
	    print STDERR "No annotation found for $line[0+$gidFlag]\n";
	}

	for (my $c = 0; $c <= $columns; $c++) {
	    print "$line[$c]" if (exists($line[$c]) && $line[$c] ne '');
	    print "\t" unless $c == $columns;
	}
	print "\n";
    }
}
