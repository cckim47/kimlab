#!/usr/bin/perl -w

use strict;

die "Usage: retrieveIDsFromEdgeCount.pl <edgeCounts.txt> <id_list>\n" unless @ARGV == 2;

open(ID,$ARGV[1]) or die "Couldn't open $ARGV[1]: $!\n";
my %idList;
while(<ID>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    $idList{$_}++;
}
close ID;

open(EDGE,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
while(<EDGE>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    next unless $idList{$line[0]};
    print join("\t",@line), "\n";
}
