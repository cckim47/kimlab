#!/usr/bin/perl -w

use strict;
use Statistics::R;

die "First run collateSTAR.pl to generate the starGenesCounts.pcl file\n" unless -e 'starstarGenesCounts.pcl';

print STDERR "Converting STAR counts\n";
system("pcl2rtable.pl starstarGenesCounts.pcl > starstarGenesCounts.txt");

my $genes = '../index/rsem/mm_grcm38/genes.gtf';

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

my $outfile = 'starstarGenesCounts_rlog.txt';
    
my $R = Statistics::R->new();
$R->startR;
$R->send(qq`library(DESeq2)`);

print STDERR "\tReading data from file\n";
$R->send(qq`data<-read.table(file="starstarGenesCounts.txt",header=T,row.names=1)`);
$R->send(qq`rounded<-round(data,0)`);
$R->send(qq`group<-c(rep("1",ceiling(ncol(data)/2)),rep("2",floor(ncol(data)/2)))`);
#$R->send(qq`group<-c(rep("1",round(ncol(data)/2),0),rep("2",ncol(data)-round(ncol(data)/2,0)))`);
$R->send(qq`frame<-data.frame(group)`);

print STDERR "\tLoading data into DESeq2\n";
$R->send(qq`dds<-DESeqDataSetFromMatrix(countData = rounded, colData = frame, design = ~ group)`);
print STDERR "\tRunning normalization\n";
$R->send(qq`rld<-rlogTransformation(dds)`);
$R->send(qq`rlogMat<-assay(rld)`);
$R->send(qq`colnames(rlogMat)<-names(data)`);
print STDERR "\tGenerating output\n";
$R->send(qq`write.table(rlogMat,file="$outfile",sep="\t",quote=FALSE,col.names=NA)`);


# add annotation and generate adj p < 0.05 file
open(DESEQ,$outfile) or die "Couldn't open $outfile for reading: $!\n";;
my $withAnn = $outfile;
$withAnn =~ s/\.txt$/_symbols\.txt/;
open(DESEQANN,">$withAnn") or die "Couldn't open $withAnn for writing: $!\n";

chomp($_=<DESEQ>); # header
my @line = split(/\t/);
$line[0] = 'UNIQID';
splice @line,1,0,'Symbol';
print DESEQANN join("\t",@line), "\n";

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
}
close DESEQ;
close DESEQANN;
