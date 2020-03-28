#!/usr/local/bin/perl
# Script to apply hapmap to a directory-full of phenotype files.

# FIXME: figure out how to get the results text file.

$|=1;				# auto-flush

use strict;
no strict "refs";

use Data::Dumper;
use Carp;

use File::Copy;   # for move 

if (scalar(@ARGV) != 1) {
    print "Usage: perl doPhenotypes.pl <max to do>\n";
}

my $maxToDo = $ARGV[0];

my @hapFiles = </DATA/susant15/input_files_withLessSubs/*.txt>;
    
foreach my $hapFile (@hapFiles) {

# get file prefix.
    $hapFile =~ /.*[\/](.*)[.][^.]*/;
    my $prefix = $1;

    my $htmlFileName = "HTML/" . $prefix . ".html";
    
    if ($maxToDo > 0) {
	my $cmd = "perl thaplomap.pl batch_query=$hapFile > $htmlFileName\n";
	print("Executing: $cmd");
	system("time $cmd");
	
	if (-e $htmlFileName) {
	    print "Probably succeeded.\n";
	    $maxToDo--;
	    if (move($hapFile, "DONE/$prefix" . ".txt") == 0) {
		print "Move to DONE failed\n.";
	    }
	}
	else {
	    print "Failed (no HTML output).\n";
	}
    }
    else {
	print "Reached maxToDo = $ARGV[0] limit.  Stopping.\n";
	exit(0);
    }
}




