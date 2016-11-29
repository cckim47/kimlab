#!/usr/bin/perl -w

use strict;

## Script to calculate performance of all workflows
## Prints performance to STDOUT

my @unitTypes = ('countsGn','countsTx','tpmGn','tpmTx','fpkmGn','fpkmTx');
my @references = ('frankenberger','haniffa','ingersoll','wong');

## Gather total number of genes assessed in each reference dataset
my %refTotals;
foreach my $ref (@references) {
    open(WCREFALL, "cat references/int_ensembl_$ref\_all.txt | wc -l |") or die;
    chomp($refTotals{$ref} = <WCREFALL>);
    close WCREFALL;
}

## Gather total number of genes called significant in each reference dataset
my %refSigCounts;
foreach my $ref (@references) {
    open(WCREFSIG, "cat references/int_ensembl_$ref\_sig.txt | wc -l |") or die;
    chomp($refSigCounts{$ref} = <WCREFSIG>);
    close WCREFSIG;
}

## Header for outfile
print "Workflow\tUnit\tReference\tPrecision\tRecall\tSigGenesOverlapWithRef\tSigGenesOverlapWithAll\n";

## Open each significant gene list, one workflow at a time
foreach my $unitType (@unitTypes) {
    my $indir = '../de/' . $unitType;
    opendir(DIR,$indir) or die "$!\n";
    while(my $file = readdir(DIR)) {
	next unless $file =~ /(.*)_sigSymbol.txt/;
	my $workflow = $1;

	print STDERR "$unitType\t$workflow\n";
	my $precisionTotal = 0;
	my $recallTotal = 0;
	my $sigCountTotal = 0;
        my $sigAllTotal = 0;
        ## Calculate performance for each reference
	foreach my $ref (@references) {
	    ## count intersection of workflow significant genes with reference significant genes
	    open(WCINTSIG, "findIntersection.pl references/int_ensembl_$ref\_sig.txt $indir/$file | sort | uniq -u | wc -l |") or die;
	    chomp(my $workflowSigIntersectionWithReferenceSig = <WCINTSIG>);
	    close WCINTSIG;
	    
	    ## count intersection of workflow significant genes with reference all represented genes
	    open(WCINTALL, "findIntersection.pl references/int_ensembl_$ref\_all.txt $indir/$file | sort | uniq -u | wc -l |") or die;
	    chomp(my $workflowSigIntersectionWithReferenceAll = <WCINTALL>);
	    close WCINTALL;

	    if ($workflowSigIntersectionWithReferenceAll == 0) {
		print STDERR "WARN: SKIPPING $workflow $unitType due to zero intersection\n";
		next;
	    }
	    
            ## Calculate performance metrics - precision and recall
	    my $precision = $workflowSigIntersectionWithReferenceSig / $workflowSigIntersectionWithReferenceAll;
	    my $recall = $workflowSigIntersectionWithReferenceSig / $refSigCounts{$ref};
	    
            ## Print results to output table
	    print "$workflow\t$unitType\t$ref";
	    print "\t$precision\t$recall\t$workflowSigIntersectionWithReferenceSig\t$workflowSigIntersectionWithReferenceAll\n";
	    
	    $precisionTotal += $precision;
	    $recallTotal += $recall;
	    $sigCountTotal += $workflowSigIntersectionWithReferenceSig;
            $sigAllTotal += $workflowSigIntersectionWithReferenceAll;
	}

        ## Calculate average performance across all four references
	next if $precisionTotal == 0;
	print "$workflow\t$unitType\taverage\t";
	print $precisionTotal/4, "\t";
	print $recallTotal/4, "\t";
	print $sigCountTotal/4, "\t";
        print $sigAllTotal/4, "\n";
    }
}
