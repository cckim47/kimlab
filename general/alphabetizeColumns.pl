#!/usr/bin/perl -w

use strict;

die "Usage: alphabetizeColumns.pl <input_pcl>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die "$!\n";

chomp(my $head = <IN>);
$head =~ s/\r//;
my @headline = split(/\t/,$head);
#print "$head\n";
shift @headline; # Oligo ID
shift @headline; # NAME
shift @headline; # GWEIGHT

# determine order of @headline data columns
my %headhash;
for my $c (0 .. $#headline) {
    $headhash{$headline[$c]} = $c;
}

# determine sort order
my @sortOrder;
foreach (sort keys %headhash) {
    push @sortOrder, $headhash{$_};
#    print "$_\t$headhash{$_}\n";
}

#print join (' ', @sortOrder);

print "UNIQID\tNAME\tGWEIGHT";
for (my $c = 0; $c <= $#sortOrder; $c++) {
    print "\t$headline[$sortOrder[$c]]";
}
print "\n";

while(<IN>) {
    chomp;
    $_ =~ s/\r//;
    if (/EWEIGHT/) {
        print "$_\n";
        next;
    }
    my @line = split(/\t/);
    print "$line[0]\t$line[1]\t1";
    for my $d ( 0 .. $#sortOrder ) {
        print "\t";
        print "$line[$sortOrder[$d]+3]" if (defined($line[$sortOrder[$d]+3]));
    }
    print "\n";
}
