#!/usr/bin/perl -w

use strict;
use Statistics::Descriptive;

# version 2

die "Usage: centerArrays.pl <input_pcl> <0 for mean, 1 for median>\n" unless @ARGV == 2;

open(IN,$ARGV[0]) or die "$!\n";

my %data; # key by column

chomp($_ = <IN>); # header
my @headSplit = split(/\t/); # headSplit
chomp($_ = <IN>); # eweight

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    my @line = split(/\t/);

    for my $c (3 .. $#line) {
	push @{$data{$c}}, $line[$c] unless $line[$c] eq '';
    }
}

my %centerValues;
foreach my $column (keys %data) {
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@{$data{$column}});
    my $center = $stat->mean if $ARGV[1] == 0;
    $center = $stat->median if $ARGV[1] == 1;
    $centerValues{$column} = $center;
}

seek IN,0,0;
chomp($_ = <IN>); # header
$_ =~ s/\r//;
print "$_\n";
chomp($_ = <IN>); # eweight
$_ =~ s/\r//;
print "$_\n";

while(<IN>) {
    chomp;
    $_ =~ s/\r//;

    my @line = split(/\t/);

    for my $c (3 .. $#line) {
	if ($line[$c] eq '') {
	    next;
	}
	else {
	    my $newVal = $line[$c] - $centerValues{$c};
	    $line[$c] = $newVal;
	}
    }

    for (my $d = 0; $d <= $#headSplit; $d++) {
	print "$line[$d]" if (exists($line[$d]) && $line[$d] ne '');
	print "\t" unless ($d == $#headSplit);
    }
    print "\n";
}
