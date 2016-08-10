#!/usr/bin/perl -w

use strict;

die "Usage: pcl2rtable.pl <pcl>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";

chomp(my $header=<IN>);
$header =~ s/\r//;
my @headSplit = split(/\t/,$header);
my $columns = $#headSplit;
print "$headSplit[0]";
for (my $c = 3; $c <= $columns; $c++) {
#    $headSplit[$c] =~ /(\d\d\d\d.b\d)/;
#    print "\t$1";
    print "\t$headSplit[$c]";
}
print "\n";

#print "GROUP";
#for (my $c = 3; $c <= $columns; $c++) {
#    print "\t1";
#}
#print "\n";

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    next if /EWEIGHT/;
    my @line = split(/\t/);
    print "$line[0]";
    for (my $c = 3; $c <= $columns; $c++) {
	if (!defined($line[$c]) || $line[$c] eq '') {
	    print "\tNA";
	}
	else {
	    print "\t$line[$c]";
	}
    }
    print "\n";
}
