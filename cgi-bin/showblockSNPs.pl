#!/usr/bin/perl
$|=1;				# auto-flush

# Given a blockSNPs file and a block number, display the SNPs in the block.

use strict;
no strict "refs";

use Data::Dumper;
use CGI ':all';
use CGI::Carp qw (warningsToBrowser fatalsToBrowser );

my %AA_CLASSES = (A => 4, B => -1, C => 3, D => 1, E => 1, F => 4, G => 2, H => 0, I => 4, J => -1, K => 0, L => 4, M => 4, N => 2, O => -1, P => 4, Q => 2, R => 0, S => 2, T => 2, U => -1, V => 4, W => 4, X => 5, Y => 2, Z => -1);

my $q = new CGI;

# When used in the hapmapper, this will be TMPDIR
# when used as part of SNPview, it will be the data directory for the data set.
my $dataPath = $q->param('dataPath');
if ($dataPath eq "") { $dataPath = 'TMPDATA'; }

my $prefix = $q->param('prefix');
my $uniquePrefix = $q->param('unique_prefix');
my $file = $q->param('file');
my $blkIdx = int($q->param('blkIdx')); # index of block (w/in chromosome)
my $blkPat = $q->param('blkPat'); # haplotype structure of block.
my $chrName = $q->param('chrName'); # chromosome
my $firstpos = int($q->param('firstpos')); # pos of first SNP in block.
my $lastpos = int($q->param('lastpos')); # pos of last SNP in block.

# collect the right SNPs out of SNPs file.
my @relevantSNPRecs = ();	# List of refs to lists
my $blockSNPsFile = "$dataPath/" . $prefix . "_SNPs.txt_chr$chrName";

# start web page here so error messages show up on it (maybe).
print header;
print start_html (-title => "SNPs in chromosome $chrName, positions $firstpos to $lastpos");

print "<STYLE TYPE=\"text/css\">\n<!--\n    TD{font-family: Arial; font-size: 10pt;\n}\n--->\n</STYLE>\n";


my $maxGenes = 0;		# maximum number of genes associated with a SNP (for table formatting).
open(SNPSFILE, "<$blockSNPsFile") || die "<pre>Open of file $blockSNPsFile failed: $!</pre>\n";
while (my $line = <SNPSFILE>) {
    chomp($line);
    my @fields = split(/\t/, $line);
    my $snpChr = $fields[0];
    my $snpPos = $fields[1] = int($fields[1]);
    my $maxFields = 0;		# max # of fields, used for number of genes.
    if (($snpChr eq $chrName) && ($snpPos >= $firstpos) && $snpPos <= $lastpos) {
	if (scalar(@fields) > $maxFields) {
	    $maxFields = scalar(@fields);
	}
	push(@relevantSNPRecs, [@fields]);
    }
    my $maxGenes = $maxFields - 4; # used for table formatting.
}

# Dumper(@relevantSNPRecs);

# Get strains, etc. from results file
my $resultsFile = "$dataPath/" . $uniquePrefix . "_results.txt";

my $phenFile = "$dataPath/" . $uniquePrefix . "_phenotypes.txt";

open(RESULTS, "<$resultsFile") || die "<pre>Internal error: results file $resultsFile open failed: $! </pre>\n";

my $datasetName = <RESULTS>;
chomp($datasetName);

my $sline = <RESULTS>;	# strain abbreviations
chomp($sline);
my @resultsStrains = split(/\t/, $sline);

close(RESULTS);

# count sizes of haplotype blocks so we can sort in decreasing order.
my @blkListPat = split('', $blkPat); # explode pattern string into list.
my %hapCounts = ();
foreach my $hap (@blkListPat) {
    if ($hap ne '?') {
	$hapCounts{$hap}++;
    }
}

# comparison function to sort INDICES of haplotype classes.
# put larger classes first, put all of same haplotype together.
sub hapcmp {
    my $hap1 = $blkListPat[$a];
    my $hap2 = $blkListPat[$b];
    my $c = ($hapCounts{$hap1} <=> $hapCounts{$hap2});
    if ($c == 0) {
	return ($hap1 cmp $hap2);
    }
    else {
	return -$c;
    }
}

# Sort strains to put haplotype blocks together.  
# CAREFUL: Strains are in different orders in the results/blocks file
# pattern and in the SNP patterns (same order as phenotype file).

# permutation of indices
my @order = sort{ hapcmp } 0 .. $#blkListPat;

# FIXME: later, put largest blocks first.
my @sortedBlkListPat = @blkListPat[@order];

# Header info
print "<h1>SNPs in chromosome $chrName, block $blkIdx, positions $firstpos to $lastpos</h1>";

# find block boundaries.
my @haploBoundaries = (0);
for (my $i= 1; $i < scalar(@sortedBlkListPat); $i++) {
    push(@haploBoundaries, (($sortedBlkListPat[$i-1] != $sortedBlkListPat[$i]) ? 1 : 0));
}

# reorder resultsStrains.
my @sortedResultsStrains = @resultsStrains[@order];

# read strain order from phenotype file (this is the order used in SNPs)
# and figure out how they need to be permuted to be in same order as
# @resultsStrains.

open(PHENFILE, "<$phenFile") || die "<pre>Internal error: results file $phenFile open failed: $! </pre>\n";

my @phenStrains = ();
while (<PHENFILE>) {
    my $line = $_;
    chomp($line);
    my ($strainAbbrev, $value) = split(/\t/, $line);
    push @phenStrains, $strainAbbrev;
}
close(PHENFILE);

# find the index permutation for the phenotype strains.
# Use sortedResultsStrains to find the destination index of each strain.
my %sortedResultsStrainIndex = ();
for (my $i = 0; $i < scalar(@sortedResultsStrains); $i++) {
    $sortedResultsStrainIndex{$sortedResultsStrains[$i]} = $i;
}

# Compute the permutation.
my @phenStrainOrder = 0 .. $#phenStrains;
for (my $i = 0; $i < scalar(@phenStrains); $i++) {
    $phenStrainOrder[$sortedResultsStrainIndex{$phenStrains[$i]}] = $i;
}

my @sortedPhenStrains = @phenStrains[@phenStrainOrder];

print "<TABLE border=0 cellspacing=3>\n";
print "<TR align=middle>\n";
# add headings for each strain abbrev
for (my $i = 0; $i < scalar(@sortedPhenStrains); $i++) {
    my $strain = $sortedPhenStrains[$i];
    my $boundary = $haploBoundaries[$i];
    if ($boundary) {
	print "<TH ALIGN=CENTER WIDTH=50>&nbsp;</TH>";
    }
    print "<TH ALIGN=CENTER WIDTH=100>$strain</TH>";
}
print "<TH>SNP ID</TH><TH>Pos</TH>";
print "<TH>Gene</TH><TH>Codon</TH></TR>\n";
foreach  my $fieldsRef (@relevantSNPRecs) {
    my ($snpChr, $snpPos, $snpID, $snpPattern, @geneInfo) = @$fieldsRef;

    # print alternating genes/codons
    print "<TR>";
    # SNP pattern
    # explode, reorder, unexplode.
    my @patList = split('', $snpPattern);
    # print "patlist = ", @patList, "\n";
    my @sortedPatList = @patList[@phenStrainOrder];

    # find major allele
    my %counts = ();
    my $major = 0;
    my $color = "white";
    
    foreach my $hap (@sortedPatList) {
	$counts{$hap}++;
    }
    if ($counts{"0"} >= $counts{1}) {
	$major = "0";
    }
    else {
	$major = "1";
    }

    # print exploded pattern.
    for (my $i = 0; $i < scalar(@sortedPatList); $i++) {
	my $haplotype = $sortedPatList[$i];
	my $boundary = $haploBoundaries[$i];
	if ($boundary) {
	    print "<TD ALIGN=CENTER WIDTH=50>&nbsp;</TD>";
	}
	if ($haplotype eq "?") {
	    $color = "white";
	}
	elsif ($haplotype eq $major) {
	    $color = "blue";
	}
	else {
	    $color = "yellow";
	}
	    
	print "<TD ALIGN=CENTER WIDTH=100 BGCOLOR=$color><FONT COLOR=$color>$haplotype</FONT></TD>";
    }

    print "<TD ALIGN=CENTER >$snpID</TD>"; # SNPID
    print "<TD ALIGN=CENTER >$snpPos</TD>";

    for (my $i = 0; $i < scalar(@geneInfo); $i = $i+2) {
	my $geneName = $geneInfo[$i];
	my $codon = $geneInfo[$i+1];
	$a = $codon;
	if ($codon =~ /NON_SYNONYMOUS_CODING/ || $codon =~ /<->/ || $codon =~ /SPLICE_SITE/) {
	    my $color = "yellow";            
            if ($a =~ /SPLICE_SITE/) {
                $a = "SPLICE_SITE";      
            }
	    elsif($a =~ /<->/) {
		$a =~ s/(.)(.)(.)\/(.)/\4/g;
		$a=~ s/([0-9]*),([A-Z])<->/\2\1/g;
		$a =~ s/!/<br \/>/g;
	    
	        while ($a =~ /(.)<->(.)/g) {
		    if($AA_CLASSES{$1} ne $AA_CLASSES{$2}) {
		        $color = "coral";
		        last;
		    }
	        } 
	    }
	    print "<TD ALIGN=CENTER BGCOLOR=$color><B><I>$geneName</I></B></TD>";
	    print "<TD ALIGN=CENTER BGCOLOR=$color><B><I>$a</I></B></TD>";
	}	
	else {
	    print "<TD ALIGN=CENTER>$geneName</TD>";
	    print "<TD ALIGN=CENTER>$codon</TD>";
	}	
    }
    print "</TR>\n";
}
print "</TABLE>\n";
print p, end_html;

close(SNPSFILE);

