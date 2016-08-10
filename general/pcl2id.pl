#!/usr/bin/perl -w

use strict;

die "Usage: pcl2id.pl <pcl>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Could not open $ARGV[0]: $!\n";

chomp($_=<IN>); # header
chomp($_=<IN>); # eweight

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    print "$line[0]\n";
}
