#!/usr/bin/perl -w

use strict;
use LWP::Simple;

die "Usage: geneFunctionalizer.pl <genelist.txt> <search string in quotes>\n" unless @ARGV == 2;

my $searchString = $ARGV[1];

my @geneList;
open(IN,$ARGV[0]) or die "Couldn't open $ARGV[0]: $!\n";
while(<IN>) {
    chomp;
    next if !$_;
    push @geneList,$_;
}
close IN;

print "Query: $searchString\n";

my %counts;
my %results;
my %hotlink;
$searchString =~ s/ /\+/g;

for my $gene (@geneList) {
#    my $query = "$gene\[Title\/Abstract\]\+AND\+\($searchString\)";
    my $query = "$gene\+AND\+\($searchString\)";

    $hotlink{$gene} = 'http://www.ncbi.nlm.nih.gov/pubmed/?term=' . $query;
    
    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url = $base . "esearch.fcgi?db=pubmed&term=$query&usehistory=n&rettype=count&retmode=text&retmax=9999";
    my $output = get($url);
    
    if ($output =~ /PhraseNotFound/) {
	$results{$gene} = 0;
	next;
    }
    elsif ($output =~ /<Count>(.*?)<\/Count>/) {
	$results{$gene} = $1;
    }
    else {
	die "Couldn't find number of results: $output\n";
    }
}

foreach my $gene (sort {$results{$b}<=>$results{$a}} keys %results) {
    print "$gene\t$results{$gene}\t$hotlink{$gene}\n";
}
