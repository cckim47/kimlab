#!/usr/bin/perl -w

use strict;

die "Usage: trimFastq5.pl <fastq> <#bases5'> <#bases3'>\n" unless @ARGV == 3;
open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
my $lineCount = 0;
while(<IN>) {
    chomp;
    next if !$_;
    next unless /^\@HWI/;
    print "$_\n";

    chomp(my $seq=<IN>);
    my $subseq = substr($seq,$ARGV[1]);
    my $len = length($subseq);
    $subseq = substr($subseq,0,$len-$ARGV[2]);
    print "$subseq\n";

    chomp($_=<IN>);
    print "$_\n";

    chomp(my $qual=<IN>);
    my $subqual = substr($qual,$ARGV[1]);
    $subqual = substr($subqual,0,$len-$ARGV[2]);
    print "$subqual\n";
}
close IN;
