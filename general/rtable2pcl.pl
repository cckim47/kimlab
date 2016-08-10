#!/usr/bin/perl -w

use strict;

die "Usage: rtable2pcl.pl <R_table_output>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
chomp(my $header = <IN>);
$header =~ s/\r//;
$header =~ s/data\.X//g;
$header =~ s/data\.//g;
$header =~ s/\" /\t/g;
$header =~ s/\"//g;
$header =~ s/\.b/\-b/g;
print "UNIQID\tNAME\tGWEIGHT\t$header\n";

print "EWEIGHT\t\t";
while($header =~ /\t/g) {
    print "\t1";
}
print "\t1\n";

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(' ');
    $line[0] =~ s/\"//g;
    print "$line[0]\t$line[0]\t1\t";
    shift @line;
    print join("\t",@line), "\n";
}    
    
