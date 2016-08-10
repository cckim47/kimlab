#!/usr/bin/perl -w

use strict;

die "Usage: trimFastq3.pl <fastq> <#bases>\n" unless @ARGV == 2;
open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
my $lineCount = 0;
while(<IN>) {
    chomp;
    next if !$_;
    next unless /^\@HWI/;
    print "$_\n";

    chomp(my $seq=<IN>);
    my $len = length($seq);
    my $subseq = substr($seq,0,$len-$ARGV[1]);
#    print "$seq\n";
    print "$subseq\n";

    chomp($_=<IN>);
    print "$_\n";

    chomp(my $qual=<IN>);
    my $subqual = substr($qual,0,$len-$ARGV[1]);
#    print "$qual\n";
    print "$subqual\n";
#    last;
}
close IN;
