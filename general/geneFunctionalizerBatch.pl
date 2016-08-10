#!/usr/bin/perl -w

use strict;

my %terms = (
    'macrophage AND polarization' => 'macPolarization',
    'macrophage AND activation' => 'macActivation',
    'diabetes OR obesity' => 'diabetesObesity',
    'atherosclerosis' => 'atherosclerosis'
    );

opendir(DIR,'.') or die "Couldn't open dir: $!\n";
while(my $file = readdir(DIR)) {
    next unless ($file =~ /(.*down)\.txt/ || $file =~ /(.*up)\.txt/);
    my $filebase = $1;

    foreach my $term (sort keys %terms) {
	print "$file against $term\n";
	my $outfile = $filebase . '_' . $terms{$term} . '.txt';
	system("./geneFunctionalizer03.pl $file \"$term\" > $outfile");
    }
}
