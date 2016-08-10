#!/usr/bin/perl -w

use strict;

die "Usage: addGenesToChart.pl <chart> <ann>\n" unless @ARGV == 2;

my %probe2symbol;
open(ANN,$ARGV[1]) or die "Couldn't open annotation: $!\n";
while(<ANN>) {
    chomp;
    $_ =~ s/\r//;
    my @line = split(/\t/);
    $probe2symbol{$line[0]} = $line[1];
}
close ANN;

open(IN,$ARGV[0]) or die "Couldn't open chart: $!\n";

chomp(my $header=<IN>);
$header =~ s/\r//;
print "$header\n";
while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    my @line = split(/\t/);
    my @ids = split(/,/,$line[5]);
    for (my $c = 0; $c <= $#ids; $c++) {
	$ids[$c] =~ s/ //g;
	if ($probe2symbol{$ids[$c]}) {
	    $ids[$c] = $probe2symbol{$ids[$c]};
	}
    }
    $line[5] = join(",",@ids);
    print join("\t",@line),"\n";
}

    
