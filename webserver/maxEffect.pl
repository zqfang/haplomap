#!/usr/bin/perl

# Script to build a table of strain counts for each file
# This is to identify examples that don't have enough strains.

# Read RESULTS directory

my @resultsFiles = <RESULTS_reprodMoreSubs/*.txt>;

foreach my $resultsFile (@resultsFiles) {

#    $resultsFile =~ /.*[\/](.*)[.][^.]*/;
#    my $prefix = $1;

    my @highEffect = ();

    open(RFILE, "<$resultsFile") || die "Open of file $resultsFile failed: $!\n";

    my $dataName = <RFILE>;	# data set
    chomp($dataName);

    my $basicName;
    my $sex;

    # remove sex suffix from data name.
    if ($dataName =~ /^(.*)_(f|m)/) {
	$basicName = $1;
	$sex = $2;
#	print "$basicName: $sex\n";
    }
    else {
	print STDERR "huh? $dataName\n";
	exit(0);
    }
	
    # read strains line
    my $strainsLine = <RFILE>;
    chomp($strainsLine);
    my @strains = split(/\t/, $strainsLine);
    my $strainCount = scalar(@strains);
#    print "$dataName\t$strainCount\n";
    <RFILE>;			# discard phenotype numbers
    
    if ($strainCount >= 12) {
	# read first 50 lines looking for effect size >= .9
	for my $i (1 .. 50) {
	    my $line =  <RFILE>;
	    chomp($line);
	    my @fields = split(/\t/, $line);
	    # Keep the line if there is a non-synonymous coding change and effect size >= .9
	    if ($fields[1] >= 0 && $fields[4] >= .85) {
		push(@highEffect, \@fields);
	    }
	}
	
	if (scalar(@highEffect) > 0) {
	    foreach my $fieldsRef (@highEffect) {
		print "$basicName\t$sex\t$fieldsRef->[0]\t$fieldsRef->[1]\t$fieldsRef->[3]\t$fieldsRef->[4]\n";
	    }
	}
    }
    close(RFILE);
}
