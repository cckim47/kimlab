#!/usr/bin/perl -w

use strict;

die "Usage: alphabetizeGenesInPCL.pl <pcl>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
chomp($_=<IN>); # header
print "$_\n";
chomp($_=<IN>); # GWEIGHT
print "$_\n";

my %data;
while(<IN>) {
    chomp;
    next if !$_;
    my @line = split(/\t/);
    $data{$line[1]} = join("\t",@line);
}
close IN;

foreach (sort keys %data) {
    print "$data{$_}\n";
}
