#!/usr/local/bin/perl


use strict;
no strict "refs";

use File::Copy;

my @resultFiles = </DATA/DillLabData/HBCGM/BATCH/HTML_lessSubs/*.html>;
my $inputFilesPath = "/DATA/susant15/input_files_withLessSubs/";

foreach my $resultFile (@resultFiles) {
	$resultFile =~ /.*[\/](.*)[.][^.]*/;
	my $filename = $1;
	if (-e $inputFilesPath.$filename.".txt") {
		print "        found it: $filename\n";
	}
	else {
  		print "moving $filename\n";
		if (move($resultFile,"DISCARD/$filename" . ".html") == 0) {
			print "MOVE of $resultFile FAILED!!!";
		}
	}
}


