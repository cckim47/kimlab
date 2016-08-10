#!/usr/bin/perl -w

use strict;

die "Usage: logTnfmTxt.pl <txt_file>)\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open Txt file: !$\n";

while(<IN>) {
  chomp;
  $_ =~ s/\r//;
  my @line = split(/\t/);
  print "$line[0]";
  for (my $c = 1; $c <= $#line; $c++) {
    print "\t";
    next if $line[$c] eq '';
    next if $line[$c] <= 0;
    print log($line[$c])/log(2);
  }
  print "\n";
}
