#!/usr/bin/perl -w

use strict;

## Script to prepare all fish (sailfish, salmon) bootstrap estimates for sleuth
## Calls an Rscript that runs wasabi for each sample, running through each sample for a given group

my @subdirs = ("SfSn4Su", "SlSl4Su", "SqSn4Su");
my @sampleBases = ("classical","nonclassical"); 
my $fragLength = 360;
my $fragSD = 200;

foreach my $subdir (@subdirs) {
	foreach my $sampleBase (@sampleBases) {
		system("Rscript model_wasabiForSu.R $subdir $sampleBase $fragLength $fragSD");
		}
	}

print("All fish samples prepped\n");
