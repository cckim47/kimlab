#!/usr/bin/perl -w

use strict;

die "Usage: meebo2entrez.pl <input_id> <annotation_table>\n" unless @ARGV == 2;

open(ANN,$ARGV[1]) or die "$!\n";

my %meebo2entrez;

while(<ANN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my ($meebo,$entrez) = split(/\t/);
#    my @temp = split(/ /,$entrez);
#    $entrez = $temp[0];
#    next if $entrez eq '-';
    $meebo2entrez{$meebo} = $entrez;
#    print "$meebo\t$entrez\n";
}

close ANN;

open(DATA,$ARGV[0]) or die "$!\n";

while(<DATA>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    if ($meebo2entrez{$_}) {
	print "$meebo2entrez{$_}\n";
    }
}
close DATA;
