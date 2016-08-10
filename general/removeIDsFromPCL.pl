#!/usr/bin/perl -w

use strict;

die "Usage: removeIDs.pl <input_pcl> <input_uniqids>" unless @ARGV == 2;

my %ids;
open(ID,$ARGV[1]) or die "$!\n";
while(<ID>) {
    chomp;
    next if !$_;
    $_ =~ s/\r//;
    $ids{$_}++;
}
close ID;

open(DATA,$ARGV[0]) or die "$!\n";
chomp(my $header = <DATA>);
print "$header\n";

while(<DATA>) {
    chomp;
    next if !$_;
    $_ =~ s/\r//;
    my @line = split(/\t/);
    next if $ids{$line[0]};
    print "$_\n";
}
close DATA;
