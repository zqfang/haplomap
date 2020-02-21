#!/usr/bin/perl
$|=1;				# auto-flush

# Form for calling SNPView.pl

use strict;
no strict "refs";

use Data::Dumper;
use CGI ':all';
use CGI::Carp qw ( fatalsToBrowser );

print header;

print start_html (-bgcolor => "#FFFAF0", 
		  -text => 'black', 
		  -link => 'blue', 
		  -vlink => 'green', 
		  -alink => 'green', 
		  -title => 'SNP viewer');

print "<h1> Combined Sanger/Peltz mouse strains SNP viewer </h1>\n";

print start_form;

print "Chromosome\t";
print popup_menu (-name => 'chromosome_name',
		  -value => ["1", "2", "3", "4", "5", "6", "7", "8", "9",
			     "10", "11", "12", "13", "14", "15", "16", "17",
			     "18", "19", "X"]);

print "Start position\t";
print textfield (-name => 'start', -size => 12);

print "End position\t";
print textfield (-name => 'end', -size => 12);

print submit(-name => 'Submit', -value => 'Submit');

print end_html;
print "\n";  # why is this needed?

