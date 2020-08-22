#!/usr/bin/perl

use strict;

# Script to build a table of strain counts for each file
# This is to identify examples that don't have enough strains.

my @resultsFiles = <RESULTS_withQuestionableSubs/*.txt>;

foreach my $resultsFile (@resultsFiles) {

#    $resultsFile =~ /.*[\/](.*)[.][^.]*/;
#    my $prefix = $1;

    my @highEffect = ();

    open(RFILE, "<$resultsFile") || die "Open of file $resultsFile failed: $!\n";

    my $dataName = <RFILE>;	# data set
    chomp($dataName);
    # read strains line
    my $strainsLine = <RFILE>;
    chomp($strainsLine);
    my @strains = split(/\t/, $strainsLine);
    my $strainCount = scalar(@strains);
    print "$dataName\t$strainCount\n";

    close(RFILE);
}
