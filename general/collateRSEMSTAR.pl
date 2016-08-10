#!/usr/bin/perl -w

use strict;

opendir(DIR,'.') or die;

my $rsemdir = '8_starrsem';
#my $gtf = '../index/rsem/mm_grcm38/genes.gtf';
my $gtf = '../index/rsem/mm_grcm38/genes.gtf';

my %geneFPKM;
my %geneTPM;
my %geneCounts;
my %geneNames;
my %isoformFPKM;
my %isoformTPM;
my %isoformCounts;
my %isoformNames;
my %samplesGene;
my %samplesIsoform;

# read GTF annotation
open(GTF,$gtf) or die "Couldn't open $gtf: $!\n";
while(<GTF>) {
    chomp;
    my @line = split(/\t/);
    if ($line[8] =~ /gene_name \"(.*?)\"/) {
	my $geneName = $1;
	$line[8] =~ /gene_id \"(.*?)\"/;
	my $geneID = $1;
	$geneNames{$geneID} = $geneName;
    }
    if ($line[8] =~ /transcript_name \"(.*?)\"/) {
        my $transcriptName = $1;
        $line[8] =~ /transcript_id \"(.*?)\"/;
        my $transcriptID = $1;
        $isoformNames{$transcriptID} = $transcriptName;
    }    
}
close GTF;

while(my $subdir = readdir(DIR)) {
    next if $subdir eq '.';
    next if $subdir eq '..';
    next unless -d $subdir;
    next unless -e "$subdir/$rsemdir";

    if (-e "$subdir/$rsemdir/$subdir.genes.results") {
	$samplesGene{$subdir}++;
	open(IN,"$subdir/$rsemdir/$subdir.genes.results") or die;
	chomp($_=<IN>); #header
	while(<IN>) {
	    chomp;
	    next if !$_;
	    my @line = split(/\t/);
	    $geneFPKM{$line[0]}{$subdir} = $line[6];
	    $geneTPM{$line[0]}{$subdir} = $line[5];
	    $geneCounts{$line[0]}{$subdir} = $line[4];
#	    $geneNames{$line[0]} = $line[4];
	}
	close IN;
    }
    
    if (-e "$subdir/$rsemdir/$subdir.isoforms.results") {
	$samplesIsoform{$subdir}++;
	open(IN,"$subdir/$rsemdir/$subdir.isoforms.results") or die;
	chomp($_=<IN>); #header
        while(<IN>) {
            chomp;
            next if !$_;
            my @line = split(/\t/);
	    my $id = "$line[0]"; 
            $isoformFPKM{$id}{$subdir} = $line[6];
            $isoformTPM{$id}{$subdir} = $line[5];
            $isoformCounts{$id}{$subdir} = $line[4];
#            $isoformNames{$id} = $line[4];
        }
        close IN;
    }


    print "$subdir\n";
}


# genes, FPKM
open(GENEPCL,'>rsemGenesFpkm.pcl') or die;
print GENEPCL "UNIQID\tGENE\tGWEIGHT";
my $sampleCount = 0;
foreach my $sampleID (sort keys %samplesGene) {
    print GENEPCL "\t$sampleID";
    $sampleCount++;
}
print GENEPCL "\n";
print GENEPCL "EWEIGHT", "\t" x 2, "\t1" x $sampleCount, "\n";

foreach my $id (sort keys %geneFPKM) {
    print GENEPCL "$id\t$geneNames{$id}\t1";
    foreach my $sampleID (sort keys %samplesGene) {
	print GENEPCL "\t$geneFPKM{$id}{$sampleID}";
    }
    print GENEPCL "\n";
}
close GENEPCL;

# genes, TPM
open(GENEPCL,'>rsemGenesTpm.pcl') or die;
print GENEPCL "UNIQID\tGENE\tGWEIGHT";
foreach my $sampleID (sort keys %samplesGene) {
    print GENEPCL "\t$sampleID";
}
print GENEPCL "\n";
print GENEPCL "EWEIGHT", "\t" x 2, "\t1" x $sampleCount, "\n";

foreach my $id (sort keys %geneFPKM) {
    print GENEPCL "$id\t$geneNames{$id}\t1";
    foreach my $sampleID (sort keys %samplesGene) {
        print GENEPCL "\t$geneTPM{$id}{$sampleID}";
    }
    print GENEPCL "\n";
}
close GENEPCL;

# genes, Counts
open(GENEPCL,'>rsemGenesCounts.pcl') or die;
print GENEPCL "UNIQID\tGENE\tGWEIGHT";
foreach my $sampleID (sort keys %samplesGene) {
    print GENEPCL "\t$sampleID";
}
print GENEPCL "\n";
print GENEPCL "EWEIGHT", "\t" x 2, "\t1" x $sampleCount, "\n";

foreach my $id (sort keys %geneFPKM) {
    print GENEPCL "$id\t$geneNames{$id}\t1";
    foreach my $sampleID (sort keys %samplesGene) {
        print GENEPCL "\t$geneCounts{$id}{$sampleID}";
    }
    print GENEPCL "\n";
}
close GENEPCL;

# isoforms, FPKM
open(ISOFORMPCL,'>rsemIsoformsFpkm.pcl') or die;
print ISOFORMPCL "UNIQID\tGENE\tGWEIGHT";
foreach my $sampleID (sort keys %samplesIsoform) {
    print ISOFORMPCL "\t$sampleID";
}
print ISOFORMPCL "\n";
print ISOFORMPCL "EWEIGHT", "\t" x 2, "\t1" x $sampleCount, "\n";

foreach my $id (sort keys %isoformFPKM) {
    print ISOFORMPCL "$id\t$isoformNames{$id}\t1";
    foreach my $sampleID (sort keys %samplesIsoform) {
        print ISOFORMPCL "\t$isoformFPKM{$id}{$sampleID}";
    }
    print ISOFORMPCL "\n";
}
close ISOFORMPCL;

# isoforms, TPM
open(ISOFORMPCL,'>rsemIsoformsTpm.pcl') or die;
print ISOFORMPCL "UNIQID\tGENE\tGWEIGHT";
foreach my $sampleID (sort keys %samplesIsoform) {
    print ISOFORMPCL "\t$sampleID";
}
print ISOFORMPCL "\n";
print ISOFORMPCL "EWEIGHT", "\t" x 2, "\t1" x $sampleCount, "\n";

foreach my $id (sort keys %isoformFPKM) {
    print ISOFORMPCL "$id\t$isoformNames{$id}\t1";
    foreach my $sampleID (sort keys %samplesIsoform) {
        print ISOFORMPCL "\t$isoformTPM{$id}{$sampleID}";
    }
    print ISOFORMPCL "\n";
}
close ISOFORMPCL;

# isoforms, Counts
open(ISOFORMPCL,'>rsemIsoformsCounts.pcl') or die;
print ISOFORMPCL "UNIQID\tGENE\tGWEIGHT";
foreach my $sampleID (sort keys %samplesIsoform) {
    print ISOFORMPCL "\t$sampleID";
}
print ISOFORMPCL "\n";
print ISOFORMPCL "EWEIGHT", "\t" x 2, "\t1" x $sampleCount, "\n";

foreach my $id (sort keys %isoformFPKM) {
    print ISOFORMPCL "$id\t$isoformNames{$id}\t1";
    foreach my $sampleID (sort keys %samplesIsoform) {
        print ISOFORMPCL "\t$isoformCounts{$id}{$sampleID}";
    }
    print ISOFORMPCL "\n";
}
close ISOFORMPCL;
