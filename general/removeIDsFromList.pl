#!/usr/bin/perl -w

use strict;

die "Usage: removeIDsFromList.pl <input_list> <IDs_to_remove>" unless @ARGV == 2;

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

while(<DATA>) {
    chomp;
    next if !$_;
    $_ =~ s/\r//;
    next if $ids{$_};
    print "$_\n";
}
close DATA;
