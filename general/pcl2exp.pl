#!/usr/bin/perl -w

use strict;

die "Usage: pcl2exp.pl <pcl>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";

chomp(my $header = <IN>);
$header =~ s/\r//;
my @headSplit = split(/\t/,$header);
shift @headSplit;
$headSplit[0] = 'AffyID';
$headSplit[1] = 'Annotation';
print join("\t",@headSplit), "\n";

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if /^EWEIGHT/;
    next if !$_;
    my @line = split(/\t/);
    splice(@line,2,1);
#    $line[0] =~ s/\W//g;
#    $line[1] =~ s/\W//g;
    $line[1] = $line[0]; # change annotation to ID to avoid failures from odd characters
    print join("\t",@line), "\n";
}
