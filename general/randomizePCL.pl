#!/usr/bin/perl -w

use strict;

die "Usage: randomizePCL.pl <pcl>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "ouldn't open $ARGV[0]: $!\n";

chomp(my $header = <IN>);
$header =~ s/\r//;
chomp(my $eweight = <IN>);
$eweight =~ s/\r//;

my @data;
while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    for (my $c = 3; $c <= $#line; $c++) {
	push @data, $line[$c];
    }
}

&fisher_yates_shuffle(\@data);

seek IN,0,0;
print "$header\n";
print "$eweight\n";
chomp($_=<IN>); # skip header
chomp($_=<IN>); # skip eweight
while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    print "$line[0]\t$line[1]\t$line[2]";
    for (my $c = 3; $c <= $#line; $c++) {
	print "\t", shift(@data);
    }
    print "\n";
}
close IN;

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
