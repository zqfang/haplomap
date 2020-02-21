#!/usr/bin/perl
$|=1;				# auto-flush

# **** To print all blocks:
# gene="*"
# pvalue = 1
# must sort by chromosome & block number
# should also format pvalue, genetic effect properly. (currently sorted by pvalue)

# Given a gene, display the haploblocks associated with the gene.

use strict;
no strict "refs";

use CGI ':all';
use CGI::Carp qw ( fatalsToBrowser );

use Digest::MD5  qw(md5 md5_hex md5_base64);

my %AA_CLASSES = (A => 4, B => -1, C => 3, D => 1, E => 1, F => 4, G => 2, H => 0, I => 4, J => -1, K => 0, L => 4, M => 4, N => 2, O => -1, P => 4, Q => 2, R => 0, S => 2, T => 2, U => -1, V => 4, W => 4, X => 5, Y => 2, Z => -1);

my $q = new CGI;

my $prefix = $q->param('prefix');
my $uniquePrefix = $q->param('unique_prefix'); 
my $set_name = $q->param("query_name");
if ($set_name eq "") {
    $set_name = "Unnamed_data_set";
}
my $data_type = $q->param("data_type");
my $p_value = $q->param("p_value");
if ($p_value eq "") {
    $p_value = 0.05;
}
my $gName = $q->param('gene_name');

my $isCategorical = ($data_type eq 'Categorical');

my $binDir = ".";	# executables for haploblocks, phmap
# my $dataPath = "/tmp";
my $dataPath = "TMPDATA";
my $haploblocksFileFull = "$dataPath/$prefix" . "_haploblocks.txt";
my $phenotypesFileFull = "$dataPath/$uniquePrefix" . "_phenotypes.txt";

# generate unique file name for gene-specific block results.
my $time = time();
my $digest = md5_hex("geneblock_" . "$time");
my $geneBlocksFileFull = "$dataPath/$digest.txt";

my @haploColors = ("red", "blue", "green", "orange", "violet", "yellow");

# copied from enhaplomap.pl
sub results_html {
    my @haploColors = ("red", "blue", "green", "orange", "violet", "yellow");
    my ($resultsFileFull, $prefix,$uniquePrefix, $p_value) = @_;
    open(RESULTS, "<$resultsFileFull") || die("<pre>Internal error: results file $resultsFileFull open failed: $! </pre>\n");

    # Dataset name
    my $datasetName = <RESULTS>;

    print header;
    print start_html (-title => "Haplotype map for gene $gName");

    chomp($datasetName);
    print "<CENTER>\n<H3>Dataset: $datasetName </H3>\n";

    # phenotype table
    print "<TABLE border=3>\n<TR align=middle>\n";
    my $sline = <RESULTS>;	# strain abbreviations
    chomp($sline);
    my @strains = split(/\t/, $sline);
    my $numStrains = scalar(@strains);
    foreach my $strain (@strains) {
	print "<TD align=\"center\" >$strain</TD>\n";
    }
    print "</TR>\n";

    # print phenotypes
    my $pline = <RESULTS>;	# phenotypes
    chomp($pline);
    my @phenotypes = split(/\t/, $pline);    
    foreach my $phenotype (@phenotypes) {
	print "<TD align=\"center\" >$phenotype</TD>\n";
    }
    print "</TR>\n";

    # print the results
    print "<H3>Significant haplotype blocks for gene $gName</H3>\n" ;
    print "<TABLE cellPadding=2 border=3>\n   <TBODY>\n  <TR align=middle>\n";
    print "<TH align=\"center\" width=\"30\">Block</TH>\n";
    if ($isCategorical) {
	print "<TH  align=\"center\" width=\"30\">F stat</TH>\n";
    }
    else {
	print "<TH  align=\"center\" width=\"30\">P-value</TH>\n";
    }
    my $hapWidth = $numStrains*3;
    print "<TH align=\"center\" width=\"30\">Genetic Effect</TH><TH  align=\"center\" colspan=$numStrains width=\"$hapWidth\">Haplotype</TH>\n";
    print "<TH align=\"center\" width=\"10\">Chr.</TH>\n<TH align=\"center\" width=\"30\">Position</TH>\n<TH align=\"center\" width=\"20\">#SNPs</TH>\n";

    # Find out maximum number of genes we're going to print.
    my $maxGenes = 0;
    while (my $bline = <RESULTS>) {
	my @fields = split(/\t/, $bline);
	my $numGenes = (scalar(@fields) - 9)/2;
	if ($numGenes > $maxGenes) {
	    $maxGenes = $numGenes;
	}
    }
    close(RESULTS);

    for (my $i = 0; $i < $maxGenes; $i++) {
	print "<TH align=\"center\" width=\"60\">Gene Name</TH>\n";
	print "<TH align=\"center\" width=\"90\">Coding?</TH>\n";
    }

    open(RESULTS, "<$resultsFileFull") || die("<pre>Internal error: results file $resultsFileFull open failed: $! </pre>\n");
    # discard the stuff we've already processed.
    <RESULTS>;			# dataset name
    <RESULTS>;			# strain abbrevs
    <RESULTS>;			# phenotype

    # now format the gene lines.
    while (my $bline = <RESULTS>) {
	chomp($bline);
	my @fields = split(/\t/, $bline);
	my $blkIdx = $fields[0];
	my $chrName = $fields[3];
	my $pattern = $fields[6];
	my $firstpos=$fields[4];
	my $lastpos=$fields[5];

	print "<TR align=middle><TD align=\"center\" ><a href=\"showblockSNPs.pl?prefix=$prefix&unique_prefix=$uniquePrefix&blkIdx=$blkIdx&blkPat=$pattern&chrName=$chrName&firstpos=$firstpos&lastpos=$lastpos\">$blkIdx</a></TD>\n";

	# p value
	print "<TD align=\"center\" >";
	printf "%.2g", $fields[7];
	print "</TD>\n";

	# effect
	print "<TD align=\"center\" >";
	printf "%.2g", $fields[8];
	print "</TD>\n";

	# print the colored haplotypes
	for (my $strIdx = 0; $strIdx < $numStrains; $strIdx++) {
	    my $allele = substr($pattern,$strIdx, 1);
	    my $color = "white";
	    if ($allele ne "?") {
		$color = $haploColors[int($allele)];
	    }
	    print "<TD align=\"center\"  BGCOLOR=$color><FONT COLOR=$color SIZE=1>$allele</TD>\n";
	}
	print "</TD>";

	print "<TD align=\"center\" >$fields[3]</TD>\n"; # chromosome name
	print "<TD align=\"center\" >$firstpos-$lastpos</TD>\n"; # position of first-last SNPs in block
	print "<TD align=\"center\" >$fields[2]</TD>\n"; # number of SNPs in block

	# starting w/ field 9, there is an alternating list of gene names/coding flags.
	# FIXME: for some reason, this is missing coding flags for genes after first in list.
	my $numGeneCols = 2*$maxGenes+9;
	for (my $i = 9; $i < $numGeneCols; $i=$i+2) {
	    if ($i < scalar(@fields)) {
		print "<TD align=\"center\" ";
		my $a = $fields[$i + 1];
		if ($fields[$i+1] ne "0") {
		    $a =~ s/(.)(.)(.)\/(.)/\4/g;
		    $a=~ s/([0-9]*),([A-Z])<->/\2\1/g;
		    $a =~ s/!/<br \/>/g;
		    my $color = "yellow";
		    while ($a =~ /(.)<->(.)/g) {
			if($AA_CLASSES{$1} ne $AA_CLASSES{$2}) {
			    $color = "coral";
			    last;
			}
		    }
		    print "BGCOLOR=$color><B><I>";
		}
		else {
		    print ">";
		}
		print "$fields[$i]"; # gene name
		if ($fields[$i+1] ne "0") {
		    print "</I></B></TD>";
		}
		else {
		    print "</TD>";
		}
		
		# coding flag
		if ($fields[$i+1] eq "1") {
		    print "<TD align=\"center\" >Y</TD>\n";
		}
		elsif($fields[$i+1] eq "0"){
		    print "<TD align=\"center\" >N</TD>\n";
		}
		else {
		    print "<TD align=\"center\" >$a</TD>\n";
		}
	    }
	    else {
		print "<TD>&nbsp;</TD><TD>&nbsp;</TD>";

	    }
	}
	print "</TR>\n";
	print "</TR>\n";
    }
    close(RESULTS);
    print p, end_html;
}


my $dtarg = "";
if ($data_type eq 'Categorical') {
    $dtarg = "-c ";
}

# build a file like the original block-oriented display, but only for the gene in question.
my $phmapcmd = "$binDir/ghmap $dtarg -l $p_value -n '$set_name' -p $phenotypesFileFull -b $haploblocksFileFull -g $gName -o $geneBlocksFileFull";

# print $phmapcmd, "\n";

system($phmapcmd);
if ($? != 0) {
    print header;

    print start_html (-title => "Haplotype map for gene $gName");
    print "<b>ghmap failed with return code $?<br>Error $!<br>";
    print "Command: $phmapcmd<br>";
    print "This is a system error -- complain to dill\@cs.stanford.edu<b><br>\n";
    print p, end_html;
    exit(1);
}

# Print the web page.

results_html($geneBlocksFileFull, $prefix, $uniquePrefix, $p_value);

# clean up temporary file.
system("rm $geneBlocksFileFull");



