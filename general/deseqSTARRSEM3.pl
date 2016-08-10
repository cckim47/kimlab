#!/usr/bin/perl -w

use strict;
use Statistics::R;

#die "First run collateSTAR.pl to generate the starGenesCounts.pcl file\n" unless -e 'starGenesCounts.pcl';

die "Usage: deseqSTARRSEM3.pl <diffexparams>\n" unless @ARGV == 1;

my $statsOutDir = 'deseqSTARRSEM3';

mkdir($statsOutDir);

#print STDERR "Converting STAR counts\n";
#system("pcl2rtable.pl starGenesCounts.pcl > deseqOnSTARcounts/merged_counts.txt");

open(PARAM,$ARGV[0]) or die;
chdir($statsOutDir);
my ($genome,$genes,$cpu);

print "Reading parameters file\n";
while(<PARAM>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    last if /^\#\#input/;
    my @line = split(/\t/);
    $genome = $line[1] if ($line[0] eq 'genome');
    $genes = $line[1] if ($line[0] eq 'genes');
    $cpu = $line[1] if ($line[0] eq 'cpu');
}


my %ann;
print STDERR "Scraping annotation from GTF\n";
open(ANN,$genes) or die "Couldn't open $genes: $!\n";
while(<ANN>) {
    chomp;
    $_ =~ s/\r//;
    my @line = split(/\t/);
    my ($id,$name);
    if (/gene_name \"(.*?)\"/) {
	$name=$1;
    }
    if (/gene_id \"(.*?)\"/) {
	$id=$1;
    }
    if ($name && $id) {
	$ann{$id} = $name;
    }
}

# read file locations
my %fileLocations;
seek PARAM,0,0;
while(<PARAM>) {
    chomp;
    last if /^\#\#inputs-starrsem/;
}
while(<PARAM>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    last if /^\#\#/;
    my @line = split(/\t/);
    $fileLocations{$line[0]} = $line[1];
}

# run comparisons
seek PARAM,0,0;
while(<PARAM>) {
    chomp;
    $_ =~ s/\r//;
    last if /^\#\#comparisons/;
}

print STDERR "Reading comparisons\n";
while(<PARAM>) {
    chomp;
    $_ =~ s/\r//;
    next if !$_;
    my @line = split(/\t/);
    my $grp1name = $line[0];
    my $grp2name = $line[1];
    my @grp1ids = split(',',$line[2]);
    my @grp2ids = split(',',$line[3]);
    my $outdir = $grp1name . 'VS' . $grp2name;
    print STDERR "Comparing $grp1name to $grp2name\n";
    next if (-e $outdir);
    mkdir($outdir);

    # create guide file and merged_counts.txt
    open(GUIDE,">$outdir/guide_file") or die "Couldn't write guide_file: $!\n";
    foreach (@grp1ids) {
        print GUIDE "../$fileLocations{$_} $_\n";
    }
    foreach (@grp2ids) {
        print GUIDE "../$fileLocations{$_} $_\n";
    }
    close GUIDE;
    system("collateRSEMusingGuideFile2.pl $outdir/guide_file > $outdir/merged_counts.txt");

    my @groupArray = ("\"$grp1name\"") x ($#grp1ids+1);
    push @groupArray, ("\"$grp2name\"") x ($#grp2ids+1);
    my $groupString = join(',',@groupArray); # NOTE: THIS ASSUMES ORDERED SAMPLES -
                                             # TRACKING OF UNORDERED SAMPLES NOT YET IMPLEMENTED
    print STDERR "\tComparison grouping: $groupString\n";

    my $sampleString = join("\",\"",@grp1ids,@grp2ids);
    $sampleString = "\"" . $sampleString . "\"";
    print "\tSamples: $sampleString\n";

    my $outfile = $grp1name . 'VS' . $grp2name . '.deseq2starrsem.results.txt';
    

    chdir($outdir);

    my $R = Statistics::R->new();
    $R->startR;
    print STDERR "\tLoading DESeq2\n";
    $R->send(qq`getwd()`);
    $R->send(qq`library(DESeq2)`);
    print STDERR "\tReading data from file\n";
    $R->send(qq`data<-read.table(file="merged_counts.txt",header=T,row.names=1)`);
    $R->send(qq`rounded<-round(data,0)`);

    $R->send(qq`write.table(rounded,file="rounded.txt",sep="\t")`);


    $R->send(qq`group<-c($groupString)`);
    $R->send(qq`samples<-c($sampleString)`);

#    $R->send(qq`group`);
#    print $R->read, "\n";

    $R->send(qq`write.table(group,file=\"group.txt\")`);

    $R->send(qq`frame<-data.frame(group,row.names=samples)`);

    $R->send(qq`write.table(frame,file=\"frame.txt\")`);

    print STDERR "\tLoading data into DESeq2\n";
    $R->send(qq`dds<-DESeqDataSetFromMatrix(countData = rounded, colData = frame, design = ~ group)`);


    print STDERR "\tRunning statistics\n";
    $R->send(qq`dds<-DESeq(dds)`);
    print STDERR "\tProducing results\n";
    $R->send(qq`res<-results(dds)`);
    print STDERR "\tOrdering by adjusted p\n";
    $R->send(qq`resOrdered<-res[order(res\$padj),]`);
    print STDERR "\tGenerating output\n";
    $R->send(qq`write.table(resOrdered,file="$outfile",sep="\t")`);
#    $R->send(qq``);

# add annotation and generate adj p < 0.05 file
    open(DESEQ,$outfile) or die "Couldn't open $outfile for reading: $!\n";;
    my $withAnn = $outfile;
    $withAnn =~ s/\.txt$/\.symbols\.txt/;
    open(DESEQANN,">$withAnn") or die "Couldn't open $withAnn for writing: $!\n";
    my $sigfile = $outfile;
    $sigfile =~ s/\.txt$/\.adjp0.05.symbols.txt/;
    open(SIG,">$sigfile") or die "Couldn't open $sigfile for writing: $!\n";
#    my $nafile = $outfile;
#    $nafile =~ s/\.txt$/\.NA.symbols.txt/;
#    open(NA,">$nafile") or die "Couldn't open $nafile for writing: $!\n";

    my $notexpressedfile = $outfile;
    $notexpressedfile =~ s/\.txt$/\.notExpressed.symbols.txt/;
    open(NOTEXPRESSED,">$notexpressedfile") or die "Couldn't open $notexpressedfile for writing: $!\n";
    
    my $outlierfile = $outfile;
    $outlierfile =~ s/\.txt$/\.outlier.symbols.txt/;
    open(OUTLIER,">$outlierfile") or die "Couldn't open $outlierfile for writing: $!\n";
    
    my $lowcountfile = $outfile;
    $lowcountfile =~ s/\.txt$/\.lowCount.symbols.txt/;
    open(LOWCOUNT,">$lowcountfile") or die "Couldn't open $lowcountfile for writing: $!\n";
    

    chomp($_=<DESEQ>); # header
    print DESEQANN "UNIQID\tSymbol\t",$_,"\n";
    print SIG "UNIQID\tSymbol\t",$_,"\n";
#    print NA "UNIQID\tSymbol\t",$_,"\n";
    print NOTEXPRESSED "UNIQID\tSymbol\t",$_,"\n";
    print OUTLIER "UNIQID\tSymbol\t",$_,"\n";
    print LOWCOUNT "UNIQID\tSymbol\t",$_,"\n";

    while(<DESEQ>) {
	chomp;
	next if !$_;
	my @line = split(/\t/);
	$line[0] =~ s/\"//g;
	if ($ann{$line[0]}) {
	    splice @line,1,0,$ann{$line[0]};
	}
	else {
	    splice @line,1,0,$line[0];
	}
	print DESEQANN join("\t",@line), "\n";

	if ($line[$#line] eq 'NA') {
#	    print NA join("\t",@line), "\n";
	    if ($line[2] == 0) {
		print NOTEXPRESSED join("\t",@line), "\n";
	    }
	    elsif ($line[5] eq 'NA') {
		print OUTLIER join("\t",@line), "\n";
	    }
	    else {
		print LOWCOUNT join("\t",@line), "\n";
	    }
	}
	elsif ($line[$#line] < 0.05) {
	    print SIG join("\t",@line), "\n";
	}
    }

    chdir('..');


}
