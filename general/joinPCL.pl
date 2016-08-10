#!/usr/bin/perl -w

#    joinPCL.pl
#
# modified from flick11.pl

use strict;

die "Usage: joinPCL.pl <file1.pcl> <file2.pcl> ...\n" unless @ARGV > 0;

&main(@ARGV);

sub main {
    my @files = @_;

    my %data = ();
    my %headers = ();
    my %names = ();
    my $colcount = 0;

    foreach (@files) {
	&readdata($_, \%data, \%names, \%headers, \$colcount);
    }

    print "UNIQID\tNAME\tGWEIGHT";
    foreach (sort {$a<=>$b} keys %headers) {
	print "\t$headers{$_}";
    }
    print "\n";

    print "EWEIGHT\t\t";
    foreach (sort {$a<=>$b} keys %headers) {
	print "\t1";
    }
    print "\n";

    foreach (keys %data) {
	my %temp = %{$data{$_}};
	print "$_\t", $names{$_}, "\t1";
#	foreach (sort {$a<=>$b} keys %temp) {
	foreach (sort {$a<=>$b} keys %headers) {
	    if (defined($temp{$_})) {
		print "\t$temp{$_}";
	    }
	    else { print "\t"; }
	}
	print "\n";
    }
}

sub readdata {
    my $file = shift;
    my $dataref = shift;
    my $nameref = shift;
    my $headerref = shift;
    my $colcountref = shift;

    open(FILE,$file) or die "Can't open PCL file\n";

    chomp(my $header=<FILE>);
    $header =~ s/\r//;
    my $gweight = 0;
    $gweight++ if $header =~ /GWEIGHT/;

    my @headercol = split(/\t/, $header);
    shift @headercol;
    shift @headercol;
    shift @headercol if $gweight == 1;

    my $startcol = $$colcountref;
    for (my $c=0; $c <= $#headercol; $c++) {
	next if !$headercol[$c];
	$headerref->{$$colcountref} = $headercol[$c];
	$$colcountref++;
    }

    while (<FILE>) {
	next if /EWEIGHT/i;
	chomp;
	$_ =~ s/\r//;
	next if !$_;
	my @line = split(/\t/);
	next if !$line[0];

	my $id = $line[0];
	$nameref->{$id} = $line[1];
	for (my $c = 2+$gweight ; $c <= $#headercol + 2 + $gweight; $c++) {
	    my $col = $startcol + $c - 2 - $gweight;
	    $dataref->{$id}{$col} = $line[$c];
	}

    }
    close FILE;
}


