#!/usr/bin/perl -wT
$|=1;				# auto-flush

# Downloads results file in spreadsheet form.
# for now, all it does is add one to the pattern digits.

use strict;
no strict "refs";

use CGI ':all';
use CGI::Carp qw ( fatalsToBrowser );

my $uniquePrefix = param("unique_prefix");
$uniquePrefix =~ /^([a-fA-F0-9]*)$/;
$uniquePrefix = $1;

my $tmpdir = "TMPDATA";
my $resultsFileFull = "$tmpdir/$uniquePrefix" . "_results.txt";

# I don't really understand this.  I wanted it to prompt for saving.
print header(-type=>'application/x-download', -attachment=>"filename=>$uniquePrefix.txt");

open(RESULTS, "<$resultsFileFull") || die "Open of $resultsFileFull failed.\n";

my $line = <RESULTS>;		# dataset name
print $line;

$line = <RESULTS>;		# strains
print $line;

$line = <RESULTS>;		# phenotype values.
print $line;

while (<RESULTS>) {
    $line = $_;
    chomp($line);
    my @fields = split(/\t/, $line);
    my @newListPat = ();
    my @listPat = unpack("C*", $fields[2]);
    foreach my $achar (@listPat) {
	if ($achar == ord('?')) {
	    push(@newListPat, $achar);
	}
	else {
	    push(@newListPat, $achar+1);
	}
    }
    my $newpat = pack("C*", @newListPat);
    $fields[2] = $newpat;
    print join("\t", @fields), "\n";
}

# print p, end_html;

close(RESULTS);
