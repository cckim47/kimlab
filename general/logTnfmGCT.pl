#!/usr/bin/perl -w

use strict;

die "Usage: logTnfmPCL.pl <gct_file>)\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open GCT file: !$\n";

chomp($_=<IN>);
$_ =~ s/\r//;
print "$_\n";

chomp($_=<IN>);
$_ =~ s/\r//;
print "$_\n";

chomp($_=<IN>);
$_ =~ s/\r//;
print "$_\n";

while(<IN>) {
  chomp;
  $_ =~ s/\r//;
  my @line = split(/\t/);
  print "$line[0]\t$line[1]";
  for (my $c = 2; $c <= $#line; $c++) {
    print "\t";
    next if $line[$c] eq '';
    next if $line[$c] <= 0;
    print log($line[$c])/log(2);
  }
  print "\n";
}
