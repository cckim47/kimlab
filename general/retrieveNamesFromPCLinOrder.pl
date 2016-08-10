#!/usr/bin/perl -w

use strict;

die "Usage: retrieveIDs.pl <input_full> <input_genelist>\nTo retrieve microarray data by ID\n" unless @ARGV == 2;

open(ID,$ARGV[1]) or die "$!\n";

my %id;
my @idOrder;
while(<ID>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    $_ =~ s/\r//;
    $id{$_}++;
    push @idOrder, $_;
}
close ID;

open(DATA,$ARGV[0]) or die "$!\n";
chomp(my $header = <DATA>);
print "$header\n";
my $tabCount = 0;
while($header =~ /\t/g) {
  $tabCount++;
}

my %data;
while(<DATA>) {
  chomp;
  $_ =~ s/\r//;
  next if !$_;
  $_ =~ s/\r//;
  if (/EWEIGHT/) {
    print "$_\n";
    next;
  }
  my $origline = $_;
  my @line = split(/\t/);
  next unless $id{$line[1]};

  my $lineTabs = 0;
  while($origline =~ /\t/g) {
    $lineTabs++;
  }
#  print "$origline";

  for (my $c = $lineTabs; $c < $tabCount; $c++) {
      $origline .= "\t";
  }
  $origline .= "\n";
  push @{$data{$line[1]}}, $origline;
}
close DATA;

for (my $c = 0; $c <= $#idOrder; $c++) {
    foreach $_ (@{$data{$idOrder[$c]}}) {
	print $_;
    }
}
