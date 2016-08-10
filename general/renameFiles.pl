#!/usr/bin/perl -w

use strict;

die "Usage: renameFiles.pl <pattern_to_change> <change_to_pattern>\n" unless @ARGV == 2;

opendir(DIR,'.') or die;

while (my $file = readdir(DIR)) {
    print "$file\n";
    next unless $file =~ /\Q$ARGV[0]\E/;
    my $newfile = $file;
    $newfile =~ s/\Q$ARGV[0]\E/$ARGV[1]/;
    print "Changing $file to $newfile\n";
    system("mv \"$file\" \"$newfile\"");
}
