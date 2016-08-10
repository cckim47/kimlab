#!/usr/bin/perl -w

use strict;

print STDERR "Usage: summarizeMapCountsSTARRSEM.pl\nRun in directory of barcoded output subdirs\n";

opendir(DIR,'.') or die "Couldn't read directory: $!\n";
my @dirs;
while(my $subdir = readdir(DIR)) {
    next if $subdir eq '.';
    next if $subdir eq '..';
    next unless -d $subdir;
    next unless -e "$subdir/8_starrsem";
    push @dirs,$subdir;
}
@dirs = sort @dirs;

my %data;
print "Sample\tMappedReads\tTotalReads\tPercMapped\n";

for (my $c = 0; $c <= $#dirs; $c++) {
    
    my $subdir = $dirs[$c];
    my $rsembam = "$subdir/8_starrsem/$subdir\.transcript.bam";
    if (-e $rsembam) {
	print STDERR "$subdir\n";
	open(STAR, "samtools view -c -F 4 $rsembam |") or die;
	my $mapped;
	while(<STAR>) {
	    chomp;
	    $mapped = $_;
	}
	close STAR;

	open(UNMAP, "samtools view -c -f 4 $rsembam |") or die;
	my $unmapped;
	while(<UNMAP>) {
	    chomp;
	    $unmapped = $_;
	}
	close UNMAP;
	print "$subdir\t$mapped\t",$unmapped+$mapped,"\t",100*$mapped/($unmapped+$mapped),"\n";
    
    }
}
