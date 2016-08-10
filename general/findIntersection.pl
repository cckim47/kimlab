#!/usr/bin/perl -w

use strict;

die "Usage: findIntersection.pl <list1> <list2> <list3> ..." unless @ARGV >= 2;

my %idTally;

foreach my $file (@ARGV) {
#    print "$file\n";
    open(FILE,$file) or die "Couldn't open $file\n";
    my %newList;
    while(<FILE>) {
        chomp;
        $_ =~ s/\r//;
        next if !$_;
	next if $newList{$_};
	$newList{$_}++;
        $idTally{$_}++;
    }
    close FILE;
}

foreach my $id (sort keys %idTally) {
    print "$id\n" if ($idTally{$id} > $#ARGV);
}
