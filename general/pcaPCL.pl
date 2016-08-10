#!/usr/bin/perl -w

use strict;
use Statistics::R;

die "Usage: pcaPCL.pl <pcl>\n" unless @ARGV == 1;

die "$ARGV[0] not found\n" unless -e $ARGV[0];

$ARGV[0] =~ /(.*)\.pcl/;
my $filebase = $1;

system("pcl2rtable.pl $ARGV[0] > $filebase\_forPCA.txt");
my $R = Statistics::R->new();
$R->send(qq`data<-read.table(file=\"$filebase\_forPCA.txt\",row.names=1,sep="\t",header=T)`);
$R->send(qq`pca<-prcomp(t(data))`);
$R->send(qq`pc<-pca\$x[,1:2]`);
$R->send(qq`png(file=\"$filebase\_forPCA.png\",width=450,height=500)`);
$R->send(qq`plot(pc[,1:ncol(pc)])`);
$R->send(qq`dev.off()`);
$R->send(qq`postscript(file=\"$filebase\_forPCA.ps\",width=5.4,height=6)`);
$R->send(qq`plot(pc[,1:ncol(pc)])`);
$R->send(qq`dev.off()`);
