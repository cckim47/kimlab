#!/usr/bin/perl -w

use strict;

die "Usage: logTnfmPCL.pl <pcl_file>)\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open PCL file: !$\n";

chomp(my $header=<IN>);
$header =~ s/\r//;
print "$header\n";

chomp(my $eweight=<IN>);
$eweight =~ s/\r//;
print "$eweight\n";

while(<IN>) {
  chomp;
  $_ =~ s/\r//;
  my @line = split(/\t/);
  print "$line[0]\t$line[1]\t$line[2]";
  for (my $c = 3; $c <= $#line; $c++) {
    print "\t";
    next if $line[$c] eq '';
    next if $line[$c] <= 0;
    print log($line[$c])/log(2);
  }
  print "\n";
}
