#!/usr/bin/perl -w

use strict;
use Statistics::Descriptive;

die "Usage: centerGenes.pl <input_pcl> <0 for mean, 1 for median>\n" unless @ARGV == 2;

open(IN,$ARGV[0]) or die "$!\n";

chomp($_ = <IN>); # header
$_ =~ s/\r//;
print "$_\n";
my @headSplit = split(/\t/);
my $columns = $#headSplit;

chomp($_ = <IN>); # eweight
$_ =~ s/\r//;
print "$_\n";

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    my @line = split(/\t/);
    my @data;
    
    for my $c (3 .. $#line) {
        if ($line[$c] eq '') {
            next;
        }
        else {
            push @data,$line[$c];
        }
    }
    
    next unless @data >= 2;
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@data);
    my $center = $stat->mean if $ARGV[1] == 0;
    $center = $stat->median if $ARGV[1] == 1;
    
    for my $c (3 .. $columns) {
        # fill in tailing tabs if missing
        unless (exists($line[$c])) {
            $line[$c] = '';
        }
        
        # skip if missing
        if ($line[$c] eq '') {
            next;
        }

        # center if numeric value present
        else {
            my $newVal = $line[$c] - $center;
            $line[$c] = $newVal;
        }
    }

    print join("\t",@line),"\n";
}
