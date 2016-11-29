#!/usr/bin/perl -w

use strict;

## Script to prepare expression matrices from expression modeling files, for use as input to differential expression tools
## Provide information about location of expression data for a given RA/EM combination in the command line call
## Expression matrix is written to STDOUT and should be redirected to a text file for saving

die "Usage: prepExpressionMatrix.pl <resultsDir> <resultsFileWithSUBDIRifNeeded> <idColumnNumberStartingFromZero> <ExpressionColumnNumberStartingFromZero> <#headerRows>\n" unless @ARGV == 5;

my ($resultsDir, $resultsFileBase, $idColumn, $expressionColumn, $headerRows) = @ARGV;
my %data;

opendir(DIR,'.') or die;

## Go into each sample directory
while(my $subdir = readdir(DIR)) {
    my $resultsFile = $resultsFileBase;
    $resultsFile =~ s/SUBDIR/$subdir/;

    next unless -e "$subdir/$resultsDir/$resultsFile";

    ## Open the file with expression data and read in data
    open(IN,"$subdir/$resultsDir/$resultsFile") or die;

    for (my $c = 1; $c <= $headerRows; $c++) {
	chomp($_=<IN>);
    }

    while(<IN>) {
	chomp;
	next unless $_;

	## Account for some lines being delimited by tabs, others by space
	my @line;
	if ($_ =~ /\t/) {
	   @line = split(/\t/);
	}
	else {
	    @line = split(/ /);
	}
	$line[$idColumn] =~ s/ //g;
	die "ID not found: $_: $!\n" unless $line[$idColumn];
	die "Expression not found in column $expressionColumn for $line[$idColumn] in $subdir/$resultsDir/$resultsFile: $!\n" unless exists($line[$expressionColumn]);

	## Clean up zero values
	$line[$expressionColumn] =~ s/^0\.000000e\+00$/0/;
	$line[$expressionColumn] =~ s/^0\.0*$/0/;
	$line[$expressionColumn] =~ s/E/e/;

	$data{$line[$idColumn]}{$subdir} = $line[$expressionColumn];
    }
    close IN;
}

## Print header for expression matrix
print "UNIQID";
foreach my $uniqid (sort keys %data) {
    foreach my $sample (sort keys %{$data{$uniqid}}) {
	print "\t$sample";
    }
    print "\n";
    last;
}

## Print expression data
foreach my $uniqid (sort keys %data) {
    print "$uniqid";
    foreach my $sample (sort keys %{$data{$uniqid}}) {
	print "\t$data{$uniqid}{$sample}";
    }
    print "\n";
}
