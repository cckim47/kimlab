#!/usr/bin/perl -w

use strict;

die "Usage: makeInformativeNames.pl <input_pcl> <parameters_file>\n" unless @ARGV == 2;

open(PARAM,$ARGV[1]) or die "$!\n";

my %names;
while(<PARAM>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
#    my $array = $line[0];
#    $names{$line[0]} = $line[1] . " \($line[0]\)";
    $names{$line[0]} = $line[1];
}

close PARAM;

open(DATA,$ARGV[0]) or die "$!\n";

chomp(my $head = <DATA>);
$head =~ s/\r//;
#print "$head\n";

my @headline = split(/\t/, $head);
my @newheadline;
foreach (@headline) {
    if ($names{$_}) {
	push @newheadline,$names{$_};
    }
    else {
	push @newheadline, $_;
    }
}

my $newhead = join("\t",@newheadline);

print "$newhead\n";

while(<DATA>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    print "$_\n";
}
