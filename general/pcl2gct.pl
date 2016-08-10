#!/usr/bin/perl -w

use strict;

die "Usage: pcl2gct.pl <gct>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Could not open $ARGV[0]: $!\n";

chomp(my $header=<IN>);
$header =~ s/\r//;
my @headSplit = split(/\t/,$header);
my $colCount = $#headSplit;

my $geneCount = 0;
while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    next if /^EWEIGHT/;
    $geneCount++;
}

print "#1.2";
for (my $c = 1; $c <= $colCount - 1; $c++) {
    print "\t";
}
print "\n";

print "$geneCount\t$colCount";
for (my $c = 2; $c <= $colCount - 1; $c++) {
    print "\t";
}
print "\n";

print "Probeset ID\tGene Symbol";
foreach (my $c = 3; $c <= $colCount; $c++) {
    print "\t$headSplit[$c]";
}
print "\n";

seek IN,0,0;
chomp($_=<IN>);
chomp($_=<IN>);
while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    print "$line[0]\t$line[1]";
    foreach (my $c = 3; $c <= $colCount; $c++) {
	print "\t$line[$c]";
    }
    print "\n";
}
