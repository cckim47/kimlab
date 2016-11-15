#!/usr/bin/perl -w

use strict;

my $indir = 'Ka';
my $infile = 'aligned.sam';
my $outdir = 'Ka4Rs';
my $outfile = 'aligned.sam';
my $conversion = '../conversion_tg.txt';


my %tid2longid;
open(ANN,$conversion) or die;
chomp($_=<ANN>); # header
while(<ANN>) {
    chomp;
    my @line = split(/\t/);
    $tid2longid{$line[0]} = $line[1] . '$$' . $line[0];
#    print "$line[0]\t$tid2longid{$line[0]}\n";
}
close ANN;


opendir(DIR,'.') or die "$!\n";
while(my $subdir = readdir(DIR)) {
    next unless -d $subdir;
    next if $subdir eq '.';
    next if $subdir eq '..';
    next if (-e "$subdir/$outdir/");
    my $read1 = $subdir . ".fastq";
    next unless (-e "$subdir/$read1");

    print "Analyzing $subdir\n";
    chdir($subdir);

    mkdir $outdir;
    open(IN,"$indir/$infile") or die;
    open(OUT,">$outdir/$outfile") or die;
    while(<IN>) {
	chomp;
	my $line = $_;
	if ($line =~ /(ENST\d+)/) {
	    my $tid = $1;
#	    print "$tid\n";
	    $line =~ s/$tid/$tid2longid{$tid}/;
	}
	print OUT "$line\n";
    }
 
    chdir('..');
}

system("notify2.pl $outdir\_MappingComplete\n");
