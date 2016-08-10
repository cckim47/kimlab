#!/usr/bin/perl -w

use strict;

die "Usage: pcl2cls.pl <pcl>\n" unless @ARGV == 1;

open(PCL,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
chomp($_=<PCL>);
$_ =~ s/\r//;
$_ =~ s/\+/p/g;
$_ =~ s/\-/n/g;

my @line = split(/\t/);
shift @line; # id
shift @line; # descript
shift @line; # gweight

my $count = 0;
my @counts = ();
my @labels = ();
my @allLabels = ();
for (my $c = 0; $c <= $#line; $c++) {
    $line[$c] =~ /(.*)\#/;
    my $base = $1;
    if ($#labels == -1) {
#	print "$c\t$count\t$base\t$allLabels[$c-1]\n";
	push @counts, $count;
	push @labels, $base;
	push @allLabels, $base;
    }
    else {
#	print "$c\t$count\t$base\t$allLabels[$c-1]\n";
	if ($base eq $allLabels[$c-1]) {
	    push @counts, $count;
	    push @allLabels, $base;
	}
	else {
	    $count++;
	    push @counts, $count;
	    push @labels, $base;
	    push @allLabels, $base;
	}
    }
}

print $#counts + 1, " ", $#labels + 1, " 1\n";
print "# ", join(" ",@labels), "\n";
print join(" ",@counts), "\n";
