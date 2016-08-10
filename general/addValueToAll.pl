#!/usr/bin/perl -w

#    goodData.pl

die "Usage: addValueToAll.pl <input_pcl> <value_to_add>\n" unless @ARGV == 2;

use strict;

open(INFILE,$ARGV[0]) or die "Can't open input file\n";

#header
chomp($_=<INFILE>);
$_ =~ s/\r//;
print "$_\n";
#eweight
chomp($_=<INFILE>);
$_ =~ s/\r//;
print "$_\n";

while (<INFILE>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    for (my $c = 3; $c <= $#line; $c++) {
	$line[$c] += $ARGV[1];
    }
    
    print join("\t",@line),"\n";
}
close INFILE;
