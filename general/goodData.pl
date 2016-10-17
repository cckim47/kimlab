#!/usr/bin/perl -w

#    goodData.pl

die "Usage: goodData.pl <input_file> <%_cutoff>\n" unless @ARGV == 2;

use strict;

my $cutoff = $ARGV[1] / 100;

open(INFILE,$ARGV[0]) or die "Can't open input file\n";

chomp(my $header=<INFILE>);
$header =~ s/\r//;
my $gweight = 0;
$gweight++ if $header =~ /GWEIGHT/i;
print "$header\n";

my @headSplit = split(/\t/,$header);
my $datacol = $#headSplit - $gweight - 1;

while (<INFILE>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;

    if (/^EWEIGHT/i) {
	print "$_\n";
	next;
    }
    
    my @line = split(/\t/);
    my $counter = 0;
    for (my $c = 2 + $gweight; $c <= $#line; $c++) {
	$counter++ if ($line[$c] ne '');
    }
    
    print "$_\n" if $counter / $datacol >= $cutoff;
}

close INFILE;
