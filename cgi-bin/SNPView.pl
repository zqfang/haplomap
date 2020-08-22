#!/usr/local/bin/perl
$|=1;				# auto-flush

# Given a blockSNPs file and a block number, display the SNPs in the block.

# This version uses sqlite

use strict;
no strict "refs";

use Data::Dumper;
use CGI ':all';
use CGI::Carp qw ( fatalsToBrowser );

use DBI;

my %AA_CLASSES = (A => 4, B => -1, C => 3, D => 1, E => 1, F => 4, G => 2, H => 0, I => 4, J => -1, K => 0, L => 4, M => 4, N => 2, O => -1, P => 4, Q => 2, R => 0, S => 2, T => 2, U => -1, V => 4, W => 4, X => 5, Y => 2, Z => -1);

my %snpHash = ();

my $q = new CGI;

# my $snpDir = $q->param('SNPdir');
my $snpDir = "SANGER_PELTZ_DATA/SNPS";
my $chrName = $q->param('chromosome_name'); # chromosome
my $firstpos = int($q->param('start')); # pos of first SNP in block.
my $lastpos = int($q->param('end')); # pos of last SNP in block.
my $rowNum = param('rownum');	# ROWID in snpAnnot table, for pagination.
my $numRowsInPage = 20;

# print "chromosome_name = $chrName\n";
# print "start = $firstpos\n";
# print "end = $lastpos\n";



# collect the right SNPs out of SNPs file.
my @relevantSNPRecs = ();	# List of refs to lists

# start web page here so error messages show up on it (maybe).
print header;

# FIXME: If database does not exist, generate an error saying it needs to be rebuilt.
my $dbh = DBI->connect("dbi:SQLite:dbname=/u01/PeltzLabData/cgi-bin/haplomap/SANGER_PELTZ_DATA/SNPGeneAnnot.db", "", "");

print "<HTML>\n<HEAD>\n<TITLE>Mouse SNP Viewer</TITLE>\n</HEAD>\n";
    
# print "<STYLE>\n";
# print "    \#vertical {\n";
# print "        width:.2em;\n";
# print "        word-wrap: break-word;\n";
# print "        align: center;\n";
# print "    }\n";
# print "</STYLE>\n";

# This was REALLY hard.
# To get vertical labels in approximately the right place, have to
# force columns to be 1em wide, and set origin at bottom and .5 em
# from the left.  (And TH & column are left aligned.)  This was
# mostly empirical -- I'm not sure if there are other workable parameters.
print "<STYLE>\n";
print "    \#vertical {\n";
print "       float:center;\n";
print "       position:relative;\n";
print "	      width:1em;\n";
print "       -webkit-transform-origin: bottom .5em;\n";
print "       -webkit-transform: rotate(-90deg);\n";
print "    }\n";
print "</STYLE>\n";


print "<title>SNP viewer: chromosome $chrName, positions $firstpos to $lastpos</title>\n";

# This is a 1 letter-wide style, used in div tags to make strain names vertical in the header table.
print "<style type=\"text/css\">\n";
print "vertical {\n";
print "width:1em;\n";
print "word-wrap: break-word;\n";
print "}\n";
print "</style>\n";

# print "<STYLE TYPE=\"text/css\">\n<!--\n    TD{font-family: Arial; font-size: 10pt;\n}\n--->\n</STYLE>\n";

# FIXME: lookup strains in table.
my @strainArrayRefs = @{ $dbh->selectall_arrayref("SELECT * FROM strains") };
my @strains = map { $_->[0] } @strainArrayRefs;
@strains = sort(@strains);

# Check whether this is a "pagination" request or an "initial" request.
if (!defined($rowNum)) {
    # it's an initial request.  Get the row number.
    # FIXME: Deal with the error if there are no rows.
    ($rowNum) = $dbh->selectrow_array("SELECT ROWID FROM snpAnnot WHERE chr = $chrName AND position >= $firstpos LIMIT 1");
}

# Retrieve the desired SNPs from the database file.

my $nextRow = $rowNum + $numRowsInPage; # first row of next page.
my $prevRow = $rowNum - $numRowsInPage; # prev row of next page.
my $endRow = $nextRow - 1;

my $snpRowsRef =
    $dbh->selectall_arrayref(
			     "SELECT ROWID, * FROM snpAnnot
				WHERE ROWID BETWEEN $rowNum AND $endRow");

# Header info
print "<h1>SNPs in chromosome $chrName</h1>\n";

# print "<TABLE align=center border=1 cellspacing=3>\n";
print "<TABLE border=1 cellspacing=3>\n";
print "<TR align=middle>\n";

print "<TH valign=\"bottom\">ID</TH><TH valign=\"bottom\">Chr</TH><TH valign=\"bottom\">Position</TH>";

# add headings for each strain abbrev
for (my $i = 0; $i < scalar(@strains); $i++) {
    my $strain = $strains[$i];
    print "<TH ALIGN=LEFT height=\"100\" valign=\"bottom\"><DIV ID=\"vertical\">$strain</DIV></TH>";
}

print "<TH valign=\"bottom\">Gene</TH><TH valign=\"bottom\">Codon</TH></TR>\n";

foreach  my $fieldsRef (@$snpRowsRef) {
    my ($rowID, $snpID, $snpChr, $snpPos, $snpAlleles, $numAnnot, $snpGene, $snpGeneAnno)
	= @$fieldsRef;

    print "<TD ALIGN=CENTER >$snpID</TD>";

    print "<TD ALIGN=CENTER >$snpChr</TD>";

    print "<TD ALIGN=CENTER >$snpPos</TD>";

    my @alleles = split('', $snpAlleles);
    foreach my $allele (@alleles) {
	print "<TD ALIGN=CENTER >$allele</TD>";
    }

    # print annotations
    print "<TD ALIGN=CENTER >$snpGene</TD>";
    print "<TD ALIGN=CENTER >$snpGeneAnno</TD>";

    print "</TR>\n";
}
print "</TABLE>\n";
print "<P>\n<H3  ALIGN=center>\n";
if ($prevRow > 0) {
    print "<A HREF=\"SNPView.pl?rownum=$prevRow\">prev $numRowsInPage</A>\n";
    # non-breaking spaces to separate prev/next.
    print "&#160;&#160;&#160;&#160;\n";
}
print "<A HREF=\"SNPView.pl?rownum=$nextRow\">next $numRowsInPage</A>\n";
print "</H3>\n";
print "</P>\n";


print p, end_html;

$dbh->disconnect();
