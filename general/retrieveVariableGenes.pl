#!/usr/bin/perl -w

use strict;
use Statistics::Descriptive;

die "Usage: retrieveMostVariableGenes.pl <txt> <percentCutoff>\n" unless @ARGV == 2;

my $total = 0;
my $input = $ARGV[0];
my $threshold = $ARGV[1];

open(IN,$input) or die "Couldn't open $input: $!\n";
chomp(my $header=<IN>); 

my %variances;
my %originalLine;
while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    my $originalLine = $_;
    my @line = split(/\t/);
    $originalLine{$line[0]} = $originalLine;
    my @data = @line[2 .. $#line];
#    print @data,"\n";
    my $stat = Statistics::Descriptive::Sparse->new;
    $stat->add_data(@data);
    my $var = $stat->variance;
    $variances{$line[0]} = $var;
    $total++;
}
close IN;

my $count = 0;
my $cutoff = int($total * $threshold / 100);
print "$header\n";
foreach my $id (sort {$variances{$b}<=>$variances{$a}} keys %variances) {
    last if $count > $cutoff;
    print "$originalLine{$id}\n";
    $count++;
}
