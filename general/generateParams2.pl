#!/usr/bin/perl -w

use strict;

die "Usage: generateParams.pl <dir1> <dir2> <dir3> ...\n" unless @ARGV >= 1;

my %fileLocations;
my %subdirPrefixes;
foreach my $dir (@ARGV) {
    opendir(DIR,$dir) or die "Couldn't open $dir for reading: $!\n";
#print $dir;
while(my $subdir = readdir(DIR)) {
	next unless -d "$dir/$subdir";
	next if $subdir eq '.';
	next if $subdir eq '..';

	opendir(SUBDIR,"$dir/$subdir") or die "Couldn't open $dir/$subdir: $!\n";
	while(my $subsubdir = readdir(SUBDIR)) {
	    if ($subsubdir =~ /1_th/) {
		die "Crash! $subdir\n" if ($fileLocations{'th'}{$subdir});
		$fileLocations{'th'}{$subdir} = $dir . '/' . $subdir . '/' . $subsubdir . '/accepted_hits.bam';
		my $prefix = $subdir;
		$prefix =~ s/\d*$//g;
		$subdirPrefixes{$prefix}{$subdir}++;
#		print "$prefix\n";
	    }
	    elsif ($subsubdir =~ /4_htseq/) {
		die "Crash! $subdir\n" if ($fileLocations{'htseq'}{$subdir});
		$fileLocations{'htseq'}{$subdir} = $dir . '/' . $subdir . '/' . $subsubdir . '/accepted_hits.counts';
		my $prefix = $subdir;
		$prefix =~ s/\d*$//g;
		$subdirPrefixes{$prefix}{$subdir}++;
#		print "$prefix\n";
	    }
	    elsif ($subsubdir =~ /5_rsem/) {
		die "Crash! $subdir\n" if ($fileLocations{'rsem'}{$subdir});
		$fileLocations{'rsem'}{$subdir} = $dir . '/' . $subdir . '/' . $subsubdir . '/' . $subdir . '.genes.results';
		my $prefix = $subdir;
		$prefix =~ s/\d*$//g;
		$subdirPrefixes{$prefix}{$subdir}++;
#		print "$prefix\n";
	    }
	    elsif ($subsubdir =~ /6_star/) {
		die "Crash! $subdir\n" if ($fileLocations{'starstar'}{$subdir});
		$fileLocations{'starstar'}{$subdir} = $dir . '/' . $subdir . '/' . $subsubdir . '/ReadsPerGene.out.tab';
		my $prefix = $subdir;
		$prefix =~ s/\d*$//g;
		$subdirPrefixes{$prefix}{$subdir}++;
#		print "$prefix\n";
	    }
	    elsif ($subsubdir =~ /7_starhtseq/) {
		die "Crash! $subdir\n" if ($fileLocations{'starhtseq'}{$subdir});
		$fileLocations{'starhtseq'}{$subdir} = $dir . '/' . $subdir . '/' . $subsubdir . '/Aligned.out.counts';
		my $prefix = $subdir;
		$prefix =~ s/\d*$//g;
		$subdirPrefixes{$prefix}{$subdir}++;
#		print "$prefix\n";
	    }
	    elsif ($subsubdir =~ /8_starrsem/) {
		die "Crash! $subdir\n" if ($fileLocations{'starrsem'}{$subdir});
		$fileLocations{'starrsem'}{$subdir} = $dir . '/' . $subdir . '/' . $subsubdir . '/' . $subdir . '.genes.results';
		my $prefix = $subdir;
		$prefix =~ s/\d*$//g;
#		print "$prefix\n";
		$subdirPrefixes{$prefix}{$subdir}++;
	    }

	}
    }
}


print "##param \n";
print "genome\t/opt/index/hs/genome.fa\n";
print "genes\t/opt/index/hs/genes.gtf\n";
print "cpu\t4\n\n";

print "##inputs-cuffdiff\n";
foreach my $sample (sort keys %{$fileLocations{'th'}}) {
    print "$sample\t$fileLocations{'th'}{$sample}\n";
}

print "\n##inputs-deseq\n";
foreach my $sample (sort keys %{$fileLocations{'htseq'}}) {
    print "$sample\t$fileLocations{'htseq'}{$sample}\n";
}

print "\n##inputs-ebseq\n";
foreach my $sample (sort keys %{$fileLocations{'rsem'}}) {
    print "$sample\t$fileLocations{'rsem'}{$sample}\n";
}

print "\n##inputs-starstar\n";
foreach my $sample (sort keys %{$fileLocations{'starstar'}}) {
    print "$sample\t$fileLocations{'starstar'}{$sample}\n";
}

print "\n##inputs-starhtseq\n";
foreach my $sample (sort keys %{$fileLocations{'starhtseq'}}) {
    print "$sample\t$fileLocations{'starhtseq'}{$sample}\n";
}

print "\n##inputs-starrsem\n";
foreach my $sample (sort keys %{$fileLocations{'starrsem'}}) {
    print "$sample\t$fileLocations{'starrsem'}{$sample}\n";
}


print "\n##comparisons\n";

#print "grp1\tgrp2\t";
#print join(',', sort keys %{$fileLocations{'starrsem'}});
#print "\t";
#print join(',', sort keys %{$fileLocations{'starrsem'}});
#print "\n";


my %comparisonTracker;
foreach my $group1 (sort keys %subdirPrefixes) {
    foreach my $group2 (sort keys %subdirPrefixes) {
	next if $group1 eq $group2;
	next if $comparisonTracker{$group1}{$group2};
	next if $comparisonTracker{$group2}{$group1};

	my @group1members = keys %{$subdirPrefixes{$group1}};
	my @group2members = keys %{$subdirPrefixes{$group2}};

	print "$group1\t$group2\t";
	print join(",", @group1members);
	print "\t";
	print join(",", @group2members);
	print "\n";

	$comparisonTracker{$group1}{$group2}++;
	$comparisonTracker{$group2}{$group1}++;
    }
}
