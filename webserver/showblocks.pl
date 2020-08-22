#!/usr/bin/perl
$|=1;				# auto-flush

# This script shows all haploblocks
# Given a gene, display the haploblocks associated with the gene.

use strict;
no strict "refs";

use CGI ':all';

use Digest::MD5  qw(md5 md5_hex md5_base64);

my $q = new CGI;

my $prefix = $q->param('prefix');
my $set_name = $q->param("query_name");
my $uniquePrefix = $q->param('unique_prefix');
if ($set_name eq "") {
    $set_name = "Unnamed_data_set";
}
my $data_type = $q->param("data_type");
my $p_value = $q->param("p_value");
if ($p_value eq "") {
    $p_value = 0.05;
}

my $isCategorical = ($data_type eq 'Categorical');

my $binDir = ".";	# executables for haploblocks, phmap
# my $dataPath = "/tmp";
my $dataPath = "TMPDATA";
my $haploblocksFileFull = "$dataPath/$prefix" . "_haploblocks.txt";
my $resultsFileFull = "$dataPath/$uniquePrefix" . "_results.txt";
my $phenotypesFileFull = "$dataPath/$uniquePrefix" . "_phenotypes.txt";
my $blockSNPsFile = "$prefix" . "_SNPs.txt";
my $blockSNPsFileFull = "$prefix" . "_SNPs.txt";

# generate unique file name for gene-specific block results.
my $time = time();
my $digest = md5_hex("geneblock_" . "$time");
my $geneBlocksFileFull = "$dataPath/$digest.txt";

my @haploColors = ("red", "blue", "green", "orange", "violet", "yellow");

# read blocks file, split lines, return a ref to an array of refs of arrays (the split lines)
sub read_blocks {
    open(HAPLOBLOCKS, "<$haploblocksFileFull") || die("<pre>Internal error: results file $haploblocksFileFull open failed: $! </pre>\n");
    my @allblocks = ();
    while (<HAPLOBLOCKS>) {
	chomp($_);
	my @blockline = split(/\t/, $_);
	push @allblocks, \@blockline;
    }
    close(HAPLOBLOCKS);
    return \@allblocks;
}

sub results_html {
    my @haploColors = ("red", "blue", "green", "orange", "violet", "yellow");
    my ($haploblocksFileFull, $blockSNPsFile, $p_value) = @_;

    open(RESULTS, "<$resultsFileFull") || die("<pre>Internal error: results file $resultsFileFull open failed: $! </pre>\n");

    # Dataset name
    my $datasetName = <RESULTS>;

    print header;
    print start_html (-title => "$datasetName");

    chomp($datasetName);
    print "<CENTER>\n<H3>Dataset: $datasetName </H3>\n";

    # phenotype table
    print "<TABLE border=3>\n<TR align=middle>\n";
    my $sline = <RESULTS>;	# strain abbreviations
    chomp($sline);
    my @strains = split(/\t/, $sline);
    my $numStrains = scalar(@strains);
    foreach my $strain (@strains) {
	print "<TD>$strain</TD>\n";
    }
    print "</TR>\n";

    # print phenotypes
    my $pline = <RESULTS>;	# phenotypes
    chomp($pline);
    my @phenotypes = split(/\t/, $pline);    
    foreach my $phenotype (@phenotypes) {
	print "<TD>$phenotype</TD>\n";
    }
    print "</TR>\n";
    print "</TABLE>\n";
    close(RESULTS);

    # read and split all the blocks, return a ref to array of refs to arrays of fields.
    my $allblocksref = read_blocks();

    # print the blocks
    print "<H3>Haplotype blocks</H3>\n" ;
    print "<TABLE cellPadding=2 border=3>\n   <TBODY>\n  <TR align=middle>\n";
    print "<TH width=\"30\">Chr.</TH>\n";
    print "<TH width=\"100\">Block</TH>\n";
    my $hapWidth = $numStrains*3;
    print "<TH width=\"200\">Position</TH>\n<TH width=\"20\">#SNPs</TH>\n";
    print "<TH colspan=$numStrains width=\"$hapWidth\">Haplotype</TH>\n";
    print "<TH width=\"20\">Coding?</TH>\n<TH width=\"60\">Genes</TH>\n</TR>\n";

    foreach my $blockref (@{$allblocksref}) {
    
	my ($chromosome, $blkidx, $start, $size, $chrbeg, $chrend, $pattern, $isCoding, @genes) = @{$blockref};

	print "<TR align=middle>\n";
	print "<TD>$chromosome</TD>\n"; # chromosome name
	print "<TD><a href=\"showblockSNPs.pl?file=$blockSNPsFile&chrName=$chromosome&blkIdx=$blkidx\">$blkidx</a></TD>\n";
	print "<TD>$chrbeg-$chrend</TD>\n"; # position of first-last SNPs in block
	print "<TD>$size</TD>\n"; # number of SNPs in block

	# print the colored haplotypes
	for (my $strIdx = 0; $strIdx < $numStrains; $strIdx++) {
	    my $allele = substr($pattern,$strIdx, 1);
	    my $color = "white";
	    if ($allele ne "?") {
		$color = $haploColors[int($allele)];
	    }
	    print "<TD BGCOLOR=$color><FONT COLOR=$color SIZE=1>$allele</TD>\n";
	}
	print "</TD>";

# Addition to show gene names	
	if ($isCoding eq "1") {
	    print "<TD>Y</TD>\n";
	}
	else {
	    print "<TD>&nbsp</TD>\n";
	}
	my $genestr = join(", ", @genes);
	if ($genestr eq "") {
	    $genestr = "&nbsp;";
	}
	print "<TD align=left width=200>$genestr</TD>\n"; 
#End of gene names call.
	print "</TR>\n";
	print "</TR>\n";
    }

    print p, end_html;
}


my $dtarg = "";
if ($data_type eq 'Categorical') {
    $dtarg = "-c ";
}

# Print the web page.

results_html($haploblocksFileFull, $blockSNPsFile, $p_value);

# clean up temporary file.
rm $geneBlocksFileFull;



