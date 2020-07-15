#!/usr/bin/perl -w
$|=1;				# auto-flush

# FIXME:  
#  * Check file size (and sanity?) for upload.
#  * Check form completeness before calling eblocks.

#  0. Figure out "categorical" from data.  Also, warn about mixtures of numerical/non-numerical data?
#  1. There may be problems with concurrent programs running.  What if phenotype files
#  are exactly the same?
#  2. Carefully think through caching, so that we don't get inappropriate cache hits.

# FLY:
#  Sort menu in numerical order.

use strict;
no strict "refs";

use IO::Handle;

use CGI ':all';
use CGI::Carp qw ( fatalsToBrowser );

use URI::Escape;

use Data::Dumper;

use Digest::MD5  qw(md5 md5_hex md5_base64);

use File::Copy;

use Time::HiRes qw(gettimeofday);
use HTML::Tooltip::Javascript;
use XML::Twig;

# System-dependent directories.
my $binDir = ".";	# executables for haploblocks, phmap

my $strains;			# ref to hash of mouse strains
my $dataDirPath;	   # directory for source-specific data files.
my $longnames;			# Prefix for generating hashes for caching 
my $geneCodingFile;		# File with gene info
my $goFileFull;			# gene ontology file.
my $expressionFileFull;		# gene expression file
my $popFileFull;           # genetic relation file 

my $haplomapURL = "http://peltz-app-02.stanford.edu/cgi-bin/haplomap/";

# batch_query should be a file name.  If it is set, we get all parameters
# from that file (no "param", web forms, etc.) and generate the output files
# without html.
my $batchQuery = param("batch_query");

# restore_parameters doesn't seem to have an error code.
if (defined $batchQuery) {
    open(BATCHQUERY, "<$batchQuery") || die "Open of file $batchQuery failed: $!\n";
    restore_parameters(\*BATCHQUERY);
    # FIXME:  Check some standard parameters to see if they were restored.
    close(BATCHQUERY);
}

my $SNPdata = param("SNPdata");
if ($SNPdata eq 'NIH') {
    $strains = {
	"A/J" => "A/J",
	"A/HeJ" => "A/H",
	"AKR/J" => "AKR",
	"C3H/HeJ" => "C3H",
	"BALB/cJ" => "B/C", 
	"BALB/cByJ" => "B/B",
	"C57BL/6J" => "C57",
	"B10.D2-H2/oSnJ" => "B10", 
	"MRL/MpJ" => "MRL",
	"NZB/BlnJ" => "NZB",
	"DBA/2J" => "DBA",
	"129/Sv" => "129",
	"NZW/LaC" => "NZW", 
	"SMJ" => "SMJ",
	"LGJ" => "LGJ",
	"LPJ" => "LPJ", 
	"BUB" => "BUB", 
	"FVB" => "FVB",
    };

    $dataDirPath = "ROCHE_DATA"; # location for Roche data files.
    $geneCodingFile = "$dataDirPath/NIH_snp_gene_coding_less10K.txt"; 
    $longnames = "ROCHE_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PERLEGEN') {

    $strains = {
	"129S1/SvImJ" => "129",
	"A/J" => "A/J",
	"AKR/J" => "AKR",
	"BALB/cByJ" => "B/B",
	"BTBR" => "BTB",
	"C3H/HeJ" => "C3H",
	"C57BL/6J" => "C57",
	"DBA/2J" => "DBA",
	"FVB/NJ" => "FVB",
	"KK/HlJ" => "KK",
	"NOD/LtJ" => "NOD",
	"NZW/LacJ" => "NZW",
    };
    
    $dataDirPath = "PERLEGEN_DATA"; # location for data files.
    $geneCodingFile = "$dataDirPath/Perlegen_snp_gene_coding_less10K.txt";
    $longnames = "PERLEGEN_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'SANGER') {
    $strains = {
	"C57" => "C57",
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"A_J" => "A/J",
	"AKR" => "AKR",
	"BALB" => "BALB",
	"C3H" => "C3H",
	"C57BL" => "C57BL",
	"CBA" => "CBA",
	"DBA" => "DBA",
	"LP_J" => "LPJ",
	"NOD" => "NOD",
	"NZO" => "NZO"
	};
    
    $dataDirPath = "SANGER_DATA"; # location for data files.
    $geneCodingFile = "$dataDirPath/sanger_gene_coding.txt";
    $longnames = "SANGER_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PSMERGED') {
    $strains = {
	"C57BL/6J" => "C57",
	"129S1/SvImJ" => "129S1",
	"A/J" => "A/J",
	"AKR/J" => "AKR",
	"C3H/HeJ" => "C3H",
	"DBA/2J" => "DBA",
	"NOD/ShiLtJ" => "NOD",
	"129P2" => "129P2",
	"129S5" => "129S5",
	"BALB/cJ" => "B/C",
	"C57BL/6NJ" => "C57BL",
	"CBA/J" => "CBA",
	"LP/J" => "LPJ",
	"NZO/HiLtJ" => "NZO",
	"BALB/cByJ" => "B/B",
	"FVB/NJ" => "FVB",
	"BTBR_T+_tf/J" => "BTB",
	"NZW/LacJ" => "NZW",
	"KK/HlJ" => "KK"
    };
    
    $dataDirPath = "PSMERGED_DATA"; # location for data files.
    $geneCodingFile = "$dataDirPath/Perlegen_snp_gene_coding_less10K.txt";
    $longnames = "PSMERGED_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'SANGER_SJL') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"BALB" => "B_C",
	"C3H" => "C3H",
	"C57BL" => "C57",
	"CAST" => "CAST",
	"CBA" => "CBA",
	"DBA" => "DBA",
	"LP_J" => "LP_J",
	"NOD" => "NOD",
	"NZO" => "NZO",
	"PWK" => "PWK",
	"SPRET" => "SPRET",
	"WSB" => "WSB",
	"SJL" => "SJL",
	};
    $dataDirPath = "SANGER_SJL_DATA"; # location for data files.
    $geneCodingFile = "$dataDirPath/sanger_gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";
    $longnames = "SANGER_SJL_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'SANGER_PELTZ') {
    $strains = {
	
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"BALB" => "B_C",
	"B10" => "B10",
	"C3H" => "C3H",
	"C57BL" => "C57",
#	"CAST" => "CAST",
	"CBA" => "CBA",
	"DBA" => "DBA",
	"FVB" => "FVB",
	"LGJ" => "LGJ",
	"LP_J" => "LP_J",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NOD" => "NOD",
	"NZB" => "NZB",
	"NZO" => "NZO",
	"NZW" => "NZW",
#	"PWK" => "PWK",
	"SJL" => "SJL",
	"SMJ" => "SMJ",
#	"SPRET" => "SPRET",
#	"WSB" => "WSB",
	};
    $dataDirPath = "SANGER_PELTZ_DATA"; # location for data files.
    $geneCodingFile = "$dataDirPath/sanger_gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";
    $longnames = "SANGER_PELTZ_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PELTZ_20121212') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"A_J" => "A/J",
	"AKR" => "AKR",
	"BALB" => "B_C",
	"C3H" => "C3H",
	"C57BL6NJ" => "C57BL6NJ",
	"CBA" => "CBA",
	"DBA" => "DBA",
	"FVB" => "FVB",
	"LPJ" => "LPJ",
	"NOD" => "NOD",
	"NZO" => "NZO",
	"B10" => "B10",
	"FVB" => "FVB",
	"LGJ" => "LGJ",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NZB" => "NZB",
	"NZW" => "NZW",
	"SMJ" => "SMJ",
	"SJL" => "SJL",
	"BTBR" => "BTBR",
	"BUB" => "BUB",
	"SWR" => "SWR",
	"SMJ" => "SMJ",
	"CAST" => "CAST",
	"PWK" => "PWK",
	"SPRET" => "SPRET",
	"WSB" => "WSB",
	};
    $dataDirPath = "PELTZ_20121212"; # location for data files.
    $geneCodingFile = "$dataDirPath/gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";
    $longnames = "PELTZ_20121212_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PELTZ_20131216') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"B10" => "B10",
	"BALB" => "BALB",
	"BTBR" => "BTBR",
	"BUB" => "BUB",
	"C3H" => "C3H",
	"C57BL6NJ" => "C57BL6NJ",
	"CAST" => "CAST",
	"CBA" => "CBA",
	"CEJ" => "CEJ",
	"DBA" => "DBA",
	"DBA1J" => "DBA1J",
	"FVB" => "FVB",
	"KK" => "KK",
	"LGJ" => "LGJ",
	"LPJ" => "LPJ",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NOD" => "NOD",
	"NON" => "NON",
	"NUJ" => "NUJ",
	"NZB" => "NZB",
	"NZO" => "NZO",
	"NZW" => "NZW",
	"PJ" => "PJ",
	"PLJ" => "PLJ",
	"PWK" => "PWK",
	"RFJ" => "RFJ",
	"RHJ" => "RHJ",
	"RIIIS" => "RIIIS",
	"SJL" => "SJL",
	"SMJ" => "SMJ",
	"SMJ" => "SMJ",
	"SPRET" => "SPRET",
	"SWR" => "SWR",
	"WSB" => "WSB",
	};
    $dataDirPath = "PELTZ_20131216"; # location for data files.
    $geneCodingFile = "$dataDirPath/gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";
    $longnames = "PELTZ_20131216_EX";
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PELTZ_20180101') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"B10" => "B10",
	"BALB" => "BALB",
	"BTBR" => "BTBR",
	"BUB" => "BUB",
	"C3H" => "C3H",
	"C57BL6NJ" => "C57BL6NJ",
	"CAST" => "CAST",
	"CBA" => "CBA",
	"CEJ" => "CEJ",
	"DBA" => "DBA",
	"DBA1J" => "DBA1J",
	"FVB" => "FVB",
	"KK" => "KK",
	"LGJ" => "LGJ",
	"LPJ" => "LPJ",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NOD" => "NOD",
	"NON" => "NON",
	"NUJ" => "NUJ",
	"NZB" => "NZB",
	"NZO" => "NZO",
	"NZW" => "NZW",
	"PJ" => "PJ",
	"PLJ" => "PLJ",
	"PWK" => "PWK",
	"RFJ" => "RFJ",
	"RHJ" => "RHJ",
	"RIIIS" => "RIIIS",
	"SJL" => "SJL",
	"SMJ" => "SMJ",
	"SMJ" => "SMJ",
	"SPRET" => "SPRET",
	"SWR" => "SWR",
	"WSB" => "WSB",
	"C57BL10J" => "C57BL10J",
	"C57BRcd" => "C57BRcd",
	"C57LJ" => "C57LJ",
	"C58" => "C58",
	"ILNJ" => "ILNJ",
	"SEA" => "SEA",
	"ST" => "ST",
	};
    $dataDirPath = "PELTZ_20180101"; # location for data files.
    $geneCodingFile = "$dataDirPath/gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";#?
    $longnames = "PELTZ_20131216_EX"; #?
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PELTZ_20190301') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"B10" => "B10",
	"BALB" => "BALB",
	"BTBR" => "BTBR",
	"BUB" => "BUB",
	"C3H" => "C3H",
	"C57BL6NJ" => "C57BL6NJ",
	"CAST" => "CAST",
	"CBA" => "CBA",
	"CEJ" => "CEJ",
	"DBA" => "DBA",
	"DBA1J" => "DBA1J",
	"FVB" => "FVB",
	"KK" => "KK",
	"LGJ" => "LGJ",
	"LPJ" => "LPJ",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NOD" => "NOD",
	"NON" => "NON",
	"NUJ" => "NUJ",
	"NZB" => "NZB",
	"NZO" => "NZO",
	"NZW" => "NZW",
	"PJ" => "PJ",
	"PLJ" => "PLJ",
	"PWK" => "PWK",
	"RFJ" => "RFJ",
	"RHJ" => "RHJ",
	"RIIIS" => "RIIIS",
	"SJL" => "SJL",
	"SMJ" => "SMJ",
	"SPRET" => "SPRET",
	"SWR" => "SWR",
	"WSB" => "WSB",
	"C57BL10J" => "C57BL10J",
	"C57BRcd" => "C57BRcd",
	"C57LJ" => "C57LJ",
	"C58" => "C58",
	"ILNJ" => "ILNJ",
	"SEA" => "SEA",
	"ST" => "ST",
	"NOR" => "NOR",
	"TALLYHO" => "TALLYHO",
	"RBF" => "RBF",
	"BPL" => "BPL",
	"BPN" => "BPN",
	};
    $dataDirPath = "PELTZ_20190301"; # location for data files.
    $geneCodingFile = "$dataDirPath/gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";#?
    $longnames = "PELTZ_20190301_EX"; #?
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'PELTZ_20200429') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"B10" => "B10",
	"BALB" => "BALB",
	"BTBR" => "BTBR",
	"BUB" => "BUB",
	"C3H" => "C3H",
	"C57BL6NJ" => "C57BL6NJ",
	"CAST" => "CAST",
	"CBA" => "CBA",
	"CEJ" => "CEJ",
	"DBA" => "DBA",
	"DBA1J" => "DBA1J",
	"FVB" => "FVB",
	"KK" => "KK",
	"LGJ" => "LGJ",
	"LPJ" => "LPJ",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NOD" => "NOD",
	"NON" => "NON",
	"NUJ" => "NUJ",
	"NZB" => "NZB",
	"NZO" => "NZO",
	"NZW" => "NZW",
	"PJ" => "PJ",
	"PLJ" => "PLJ",
	"PWK" => "PWK",
	"RFJ" => "RFJ",
	"RHJ" => "RHJ",
	"RIIIS" => "RIIIS",
	"SJL" => "SJL",
	"SMJ" => "SMJ",
	"SPRET" => "SPRET",
	"SWR" => "SWR",
	"WSB" => "WSB",
	"C57BL10J" => "C57BL10J",
	"C57BRcd" => "C57BRcd",
	"C57LJ" => "C57LJ",
	"C58" => "C58",
	"ILNJ" => "ILNJ",
	"SEA" => "SEA",
	"ST" => "ST",
	"NOR" => "NOR",
	"TALLYHO" => "TALLYHO",
	"RBF" => "RBF",
	"BPL" => "BPL",
	"BPN" => "BPN",
	};
    $dataDirPath = "PELTZ_20200429"; # location for data files.
    $geneCodingFile = "$dataDirPath/gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";#?
    $longnames = "PELTZ_20200429_EX"; #?
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
	$popFileFull="$dataDirPath/mouse54_grm.rel"
}
elsif ($SNPdata eq 'PELTZ_20200505') {
    $strains = {
        "C57BL/6J" => "C57/6J", #reference
	"129P2" => "129P2",
	"129S1" => "129S1",
	"129S5" => "129S5",
	"AKR" => "AKR",
	"A_J" => "A/J",
	"B10" => "B10",
	"BALB" => "BALB",
	"BTBR" => "BTBR",
	"BUB" => "BUB",
	"C3H" => "C3H",
	"C57BL6NJ" => "C57BL6NJ",
	"CAST" => "CAST",
	"CBA" => "CBA",
	"CEJ" => "CEJ",
	"DBA" => "DBA",
	"DBA1J" => "DBA1J",
	"FVB" => "FVB",
	"KK" => "KK",
	"LGJ" => "LGJ",
	"LPJ" => "LPJ",
	"MAMy" => "MAMy",
	"MRL" => "MRL",
	"NOD" => "NOD",
	"NON" => "NON",
	"NUJ" => "NUJ",
	"NZB" => "NZB",
	"NZO" => "NZO",
	"NZW" => "NZW",
	"PJ" => "PJ",
	"PLJ" => "PLJ",
	"PWK" => "PWK",
	"RFJ" => "RFJ",
	"RHJ" => "RHJ",
	"RIIIS" => "RIIIS",
	"SJL" => "SJL",
	"SMJ" => "SMJ",
	"SPRET" => "SPRET",
	"SWR" => "SWR",
	"WSB" => "WSB",
	"C57BL10J" => "C57BL10J",
	"C57BRcd" => "C57BRcd",
	"C57LJ" => "C57LJ",
	"C58" => "C58",
	"ILNJ" => "ILNJ",
	"SEA" => "SEA",
	"ST" => "ST",
	"NOR" => "NOR",
	"TALLYHO" => "TALLYHO",
	"RBF" => "RBF",
	"BPL" => "BPL",
	"BPN" => "BPN",
	"BALBBYJ" => "BALBBYJ",
	};
    $dataDirPath = "PELTZ_20200505"; # location for data files.
    $geneCodingFile = "$dataDirPath/gene_coding.txt";
    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";#?
    $longnames = "PELTZ_20200505_EX"; #?
    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}
elsif ($SNPdata eq 'FLY') {
    $strains = {
        "REF" => "REF", #reference
	"S21" => "S21",
	"S26" => "S26",
	"S28" => "S28",
	"S38" => "S38",
	"S40" => "S40",
	"S41" => "S41",
	"S42" => "S42",
	"S45" => "S45",
	"S49" => "S49",
	"S57" => "S57",
	"S59" => "S59",
	"S69" => "S69",
	"S73" => "S73",
	"S75" => "S75",
	"S83" => "S83",
	"S85" => "S85",
	"S88" => "S88",
	"S91" => "S91",
	"S93" => "S93",
	"S101" => "S101",
	"S105" => "S105",
	"S109" => "S109",
	"S129" => "S129",
	"S136" => "S136",
	"S138" => "S138",
	"S142" => "S142",
	"S149" => "S149",
	"S153" => "S153",
	"S158" => "S158",
	"S161" => "S161",
	"S176" => "S176",
	"S177" => "S177",
	"S181" => "S181",
	"S195" => "S195",
	"S208" => "S208",
	"S217" => "S217",
	"S227" => "S227",
	"S228" => "S228",
	"S229" => "S229",
	"S233" => "S233",
	"S235" => "S235",
	"S237" => "S237",
	"S239" => "S239",
	"S256" => "S256",
	"S272" => "S272",
	"S280" => "S280",
	"S287" => "S287",
	"S309" => "S309",
	"S310" => "S310",
	"S313" => "S313",
	"S317" => "S317",
	"S318" => "S318",
	"S320" => "S320",
	"S321" => "S321",
	"S325" => "S325",
	"S332" => "S332",
	"S338" => "S338",
	"S350" => "S350",
	"S352" => "S352",
	"S356" => "S356",
	"S357" => "S357",
	"S358" => "S358",
	"S359" => "S359",
	"S362" => "S362",
	"S365" => "S365",
	"S367" => "S367",
	"S370" => "S370",
	"S371" => "S371",
	"S373" => "S373",
	"S374" => "S374",
	"S375" => "S375",
	"S377" => "S377",
	"S378" => "S378",
	"S379" => "S379",
	"S380" => "S380",
	"S381" => "S381",
	"S383" => "S383",
	"S386" => "S386",
	"S391" => "S391",
	"S392" => "S392",
	"S398" => "S398",
	"S399" => "S399",
	"S405" => "S405",
	"S406" => "S406",
	"S409" => "S409",
	"S426" => "S426",
	"S427" => "S427",
	"S437" => "S437",
	"S439" => "S439",
	"S440" => "S440",
	"S441" => "S441",
	"S443" => "S443",
	"S461" => "S461",
	"S491" => "S491",
	"S492" => "S492",
	"S502" => "S502",
	"S508" => "S508",
	"S509" => "S509",
	"S513" => "S513",
	"S517" => "S517",
	"S531" => "S531",
	"S535" => "S535",
	"S554" => "S554",
	"S555" => "S555",
	"S563" => "S563",
	"S589" => "S589",
	"S591" => "S591",
	"S595" => "S595",
	"S639" => "S639",
	"S642" => "S642",
	"S646" => "S646",
	"S703" => "S703",
	"S705" => "S705",
	"S707" => "S707",
	"S712" => "S712",
	"S714" => "S714",
	"S716" => "S716",
	"S721" => "S721",
	"S727" => "S727",
	"S730" => "S730",
	"S732" => "S732",
	"S737" => "S737",
	"S738" => "S738",
	"S757" => "S757",
	"S761" => "S761",
	"S765" => "S765",
	"S774" => "S774",
	"S776" => "S776",
	"S783" => "S783",
	"S786" => "S786",
	"S787" => "S787",
	"S790" => "S790",
	"S796" => "S796",
	"S799" => "S799",
	"S801" => "S801",
	"S802" => "S802",
	"S804" => "S804",
	"S805" => "S805",
	"S808" => "S808",
	"S810" => "S810",
	"S812" => "S812",
	"S818" => "S818",
	"S820" => "S820",
	"S822" => "S822",
	"S832" => "S832",
	"S837" => "S837",
	"S852" => "S852",
	"S855" => "S855",
	"S857" => "S857",
	"S859" => "S859",
	"S861" => "S861",
	"S879" => "S879",
	"S882" => "S882",
	"S884" => "S884",
	"S887" => "S887",
	"S890" => "S890",
	"S892" => "S892",
	"S894" => "S894",
	"S897" => "S897",
	"S907" => "S907",
	"S908" => "S908",
	"S911" => "S911",
    };
    $dataDirPath = "FLY_DATA"; # location for data files.
    $geneCodingFile = "$dataDirPath/fly_snp_coding.txt";
#    $goFileFull = "$dataDirPath/genes_to_go_terms.txt";
    $longnames = "FLY_DATA";
#    $expressionFileFull = "$dataDirPath/" . "compact_gene_expr.txt";
}


my $SNPSdataDirPath = "$dataDirPath/SNPS"; # location for chromosome SNPs data files.

# my $dataDirPath = "PERLEGEN_DATA"; # location for data files.

# where to store temporary files
my $tmpdir = "TMPDATA";
# my $tmpdir = "/tmp";

my $images_url = "http://gmc.stanford.edu/icons/images";

my $run = param ('run');


if ( $run eq 'Find' ) {
    &output_display;
}
elsif ($run eq 'Save Query') {
    # save a copy of the form
    my $set_name = param("query_name");
    if ($set_name eq "") {
	$set_name = "Unnamed_data_set";
    }
    print header(-type=>'application/x-download', -attachment=>"$set_name.txt");
#    print header, start_html;
#    print "<pre>";
    save_parameters( \*STDOUT );
#    print "</pre>";
#    print end_html;
}
elsif ($run eq 'Upload Form Fields') {
    # populate menu from uploaded file.

    print header;

    print start_html (-bgcolor => "#FFFAF0", 
		      -text => 'black', 
		      -link => 'blue', 
		      -vlink => 'green', 
		      -alink => 'green', 
		      -title => 'msnpdb Stanford Medical Center');

    my $lightweight_fh = upload('uploaded_form');

    if (!defined $lightweight_fh) {
	# complain and redisplay form.
	print "<h2>Bad upload file.  Please try again.</h2>\n";
	&print_form;
	exit(0);
    }
    else {
	# Upgrade the handle to one compatible with IO::Handle:
	restore_parameters(*$lightweight_fh);
	&print_form;
    }
}
else {

    print header;
    
    print start_html (-bgcolor => "#FFFAF0", 
		      -text => 'black', 
		      -link => 'blue', 
		      -vlink => 'green', 
		      -alink => 'green', 
		      -title => 'msnpdb Stanford Medical Center');
    
    &print_form;
}

exit(0);


# cut-and-pasted from a perl FAQ.  Make numbers into a string with commas
# separating groups of three digits.
sub commify {
    local $_ = shift;
    s{(?:(?<=^)|(?<=^-))(\d{4,})}
    {my $n = $1;
     $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
		  $n;
	      }e;
     return $_;
 }

# comparison with alphabetic prefixes, but in numerical order.
sub anumcmp {
    if (! ($a =~ /(\D*)(\d*)/)) {
	print "Why did that not match?!?\n";
	};
    my $alpha_a = $1;
    my $num_a = $2;
    $b =~ /(\D*)(\d*)/;
    my $alpha_b = $1;
    my $num_b = $2;
    my $cmp1 = ($alpha_a cmp $alpha_b);
    if ($cmp1 == 0) {
	return ($num_a <=> $num_b);
    }
    else {
	return $cmp1;
    }
}

#populates genesToMesh so that it is a mapping from gene or protein product names to sets of mesh terms
sub enterMappings {
    my ($twig,$descrip,$genesToMeshRef) = @_;
    my %meshSet = ();
    foreach my $meshTerm ($descrip->first_child('MeshTermList')->children) {
        $meshSet{$meshTerm->text} = 1;
    }
    $genesToMeshRef->{$descrip->first_child('Name')->text} = \%meshSet;
    $descrip->purge();
}

#returns the tooltip text for the given gene
sub tooltip_text {
    my $text = "";
    my ($geneName,$genesToMeshRef) = @_; 
    foreach my $meshTerm (keys (%{$genesToMeshRef->{lc($geneName)}})) {
        $text = $text." ".$meshTerm;
    }
    return($text);
}

# Uploads a bunch of param values from a file.
# This is used for a form or for batch execution of files.
# CGI's "lightweight filehandles" are very badly documented,
# so I'm guessing that the argument can either be a string
# file name, or result of "upload" from a form field.



sub results_html {
    my @haploColors = ("red", "blue", "green", "orange", "violet", "yellow");
    # my @tisabbrevs = ("BC", "CB", "CD", "HP", "KD", "LV", "LN", "PF", "QD", "SC", "SP", "SA", "TC");
    # my @tissues = ("B Cell Basal", "Cerebellum", "Chondrocyte", 
	# 	   "Hippocampus", "Kidney", "Liver", "Lung", "Prefrontal Cortex",
	# 	   "Quadricep", "Spinal Cord", "Spleen", "Striatum", "T Cell Basal");

	my @tisabbrevs = ("AOR", "ACO", "BSK", "BMA", "CBM", "CCX", "DCN", "DPM", "ENP", "EXP", "GFP", "HEA", "HLA", "HLRA", "HLV","HRA","HRV", 
	                  "HIP","IBA","KNE","LMU","LIV","LUN","MAG","MAT","SPL", "STR","SAT","THY", "TOG","TRC", "URB");

	my @tissues	= ("aorta","ascending colon","back skin","bone marrow","cerebellum","cerebral cortex","descending colon",
			"diaphragm","endocrine pancreas","exocrine pancreas","gonadal fat pad","heart","heart left atrium",
			"heart left atrium and heart right atrium","heart left ventricle","heart right atrium","heart right ventricle",
			"hippocampus","interscapular brown adipose tissue","kidney","limb muscle","liver","lung","mammary gland",
			"mesenteric adipose tissue","spleen","striatum","subcutaneous adipose tissue","thymus",
			"tongue","trachea","urinary bladder");
				
    
    my ($data_type, $resultsFileFull,  $prefix, $uniquePrefix, $p_value, @regionChoices) = @_;
    
    my $isCategorical = ($data_type eq 'Categorical');
    
    open(RESULTS, "<$resultsFileFull") || confess "Internal error: results file $resultsFileFull open failed: $!\n";
    
    # Dataset name
    my $datasetName = <RESULTS>;	# need to untaint this.
    chomp($datasetName);
    $datasetName =~ /(.*)/;
    $datasetName = $1;
    
    print "<HTML>\n<HEAD>\n<TITLE>$datasetName</TITLE>\n</HEAD>\n";
    
    print "<STYLE>\n";
    print "    \#download {\n";
    print "        float:left;\n";
    print "    }\n";
    
    print "    \#showblocks {\n";
    print "        float:right;\n";
    print "    }\n";
    
    print "    \#legend {\n";
    print "        clear:both;\n";
    print "        margin-top:50;\n";
    print "    }\n";
    print "</STYLE>\n";
    
    print "<BODY>\n";
    
#    print "<h2> SNPdata = $SNPdata </h2>\n";
    
    print "<CENTER>\n<H3>Dataset: $datasetName </H3></CENTER>\n";
    
    my $download_file = "$haplomapURL" . "downloadResults.pl?unique_prefix=$uniquePrefix";
    
#     print "<TABLE CELLSPACING=100>\n<TR>\n";	# some random commands
#     print "<TR><TD><a href=\"$download_file\">Results download</a></TD>\n";
#     print "<TD><a href=\"showblocks.pl?prefix=$prefix\">Show blocks (warning: big!)</a></TD>\n"; 
#     print "</TABLE><P>\n";
    
    print "<DIV ID=\"download\">\n";
    print "<a href=\"$download_file\">Results download</a>\n";
    print "</DIV>\n";
    
    print "<DIV ID=\"showblocks\">\n";
    my $showblocksScript = $haplomapURL . "showblocks.pl?prefix=$prefix"."&unique_prefix=$uniquePrefix";
    print "<a href=\"$showblocksScript\">Show blocks (warning: big!)</a>\n"; 
    print "</DIV>\n";
    
    # Print gene expression legend only fo Mouse genome
    if (defined $expressionFileFull) {
	print "<DIV ID=\"legend\">\n";
	print "<H3>Gene Expression Legend</H3>\n";
	
	print "<TABLE>\n<TR>\n";	# outer table for all legend stuff.
	print "<TD  VALIGN=top WIDTH=200>\n";		# first sub table -- color code.
	
	print "<TABLE WIDTH=\"100\">\n";
	print "<TR><TD>Present</TD><TD WIDTH=\"20\" BGCOLOR=blue></TD></TR>\n";
	print "<TR><TD>Absent</TD> <TD WIDTH=\"20\" BGCOLOR=green></TD></TR>\n";
	print "<TR><TD>No data</TD> <TD WIDTH=\"20\" BGCOLOR=#BEBEBE></TD></TR>\n";
	print "</TABLE>";
	
	print "<DIV ID=\"tissues\">\n";
	# next n columns -- tissue abbreviations.
	print "</TD>\n";
	
	for (my $i = 0; $i < scalar(@tisabbrevs); $i++) {
	    print "</TABLE></TD>" if ($i > 0 && $i % 4 == 0);
	    print "<TD  VALIGN=top WIDTH=200><TABLE BORDER=3 WIDTH=100>" if ($i % 4 == 0);
	    print "<TR><TD>$tisabbrevs[$i]</TD><TD>$tissues[$i]</TD></TR>\n";
	}
	print "</TABLE>\n";
	
	print "</TD></TR></TABLE>\n"; # finish outer table for legend
	print "</DIV>\n";		# end legend.
    }
    
    # phenotype table
    print "<H3>Phenotypes</H3>\n";
    print "<TABLE border=3>\n<TR align=middle>\n";
    my $sline = <RESULTS>;	# strain abbreviations
    chomp($sline);
    my @strains = split(/\t/, $sline);
    my $numStrains = scalar(@strains);
    foreach my $strain (@strains) {
	print "<TH>$strain</TH>\n";
    }
    print "</TR>\n";
    
    # print phenotypes
    my $pline = <RESULTS>;	# phenotypes
    chomp($pline);
    my @phenotypes = split(/\t/, $pline);    
    foreach my $phenotype (@phenotypes) {
	print "<TD>$phenotype</TD>\n";
    }
    print "</TR></TABLE>\n";
    
    # print the results
    print "<H3>Genes with significant matching haplotype blocks</H3>\n" ;
    print "<TABLE cellPadding=2 border=3>\n   <TBODY>\n  <TR align=middle>\n";
    print "<TH width=\"200\">Gene</TH>\n";
    if ($isCategorical) {
	print "<TH width=\"100\">F stat</TH>\n";
    }
    else {
	print "<TH width=\"100\">P-value</TH>\n";
    }
    my $hapWidth = $numStrains*7 + 2; # I think: 5 per color, one extra cell padding of 2.
    print "<TH width=\"100\">Genetic Effect</TH><TH width=\"100\">FDR</TH>\n";
	print "<TH width=\"100\">pop P-value</TH><TH width=\"100\">pop FDR</TH><TH colspan=$numStrains width=\"$hapWidth\">Haplotype</TH>\n";
    print "<TH width=\"10\">Chr.</TH>\n<TH width=\"200\">Position</TH>\n\n";

    my $numtissues = scalar(@tisabbrevs);
    
    if (defined $expressionFileFull) {
	for my $tis (@tisabbrevs) {
	    print "<TD width=10>$tis</TD>";
	}
	print "</TR>\n";
    }

    # This is for "Show only genes expressed in ..."
    my @tissueIndices = ();
    my %tissueMap = map { $_ => 1 } @regionChoices;
    
    for(my $i = 0; $i < @tissues; $i++) {
	if (exists($tissueMap{$tissues[$i]})) {
	    @tissueIndices = (@tissueIndices, $i);
	}
    }
    my %genesToMesh = ();
###    my $t = XML::Twig->new (twig_handlers =>
###                            { Descriptor => sub{ enterMappings(@_,\%genesToMesh); } }
###                        );
###    $t->parsefile("forScripts/genes.xml");
###    $t->purge();

    my $tt = HTML::Tooltip::Javascript->new(
        javascript_dir => '/u01/PeltzLabData/javascripts/',
        options => {
            default_tip => 'No mesh terms found',
            delay => 1
        },
    );

    while (my $bline = <RESULTS>) {
	my @fields = split(/\t/, $bline);
	my $present = $fields[12];
	
	my $skip = 0;
	foreach my $index (@tissueIndices){
	    $skip = $skip || substr($present, $index, 1) !~ /[P-]/;
	}
	
	next if ($skip);
	
	# print gene name, bold & italic if it's coding.
	print "<TR align=middle><TD";
	if ($fields[1] eq "0") {
	    print " BGCOLOR=yellow><B><I>";
	}
	elsif ($fields[1] eq "1") {
	    print " BGCOLOR=coral><B><I>";
	}
        elsif ($fields[1] eq "2") {
            print " BGCOLOR=LightSkyBlue><B><I>";
        }
	else {
	    print ">";
	}
	
	# watch out for spaces, etc. in data set name
	my $datasetName_safe = uri_escape($datasetName);

	my $showgeneblocksScript = $haplomapURL . "showgeneblocks.pl?prefix=$prefix&unique_prefix=$uniquePrefix&query_name=$datasetName_safe&gene_name=$fields[0]&data_type=$data_type&p_value=$p_value";
	my $meshTerms = tooltip_text($fields[0],\%genesToMesh);
	my $meshTooltip = $tt->tooltip($meshTerms);
	print "<a href=$showgeneblocksScript $meshTooltip>";
	print "$fields[0]</a>";

	if ($SNPdata eq 'FLY') {
	    # this doesn't work, unfortunately.
	    print " <a href=\"http://flybase.net/.bin/fbgenq.html?symbol=$fields[0]\">*</a>";
	}
	else {
	    print " <a href=\"http://www.informatics.jax.org/searchtool/Search.do?query=$fields[0]\">*</a>";
	}
	if ($fields[1] ne "-1") {
	    print "</B></I>";
	}
	print "</TD>\n";		
   # p-value
	print "<TD>";
	printf "%.2g", $fields[3];
	print "</TD>\n";
	
	# effect
	print "<TD>";
	printf "%.2g", $fields[4];
	print "</TD>\n";

	# FDR
	print "<TD>";
	printf "%.2g", $fields[5];
	print "</TD>\n";
	# popPvalue
	print "<TD>";
	printf "%.2g", $fields[6];
	print "</TD>\n";
	# popFDR
	print "<TD>";
	printf "%.2g", $fields[7];
	print "</TD>\n";
	# print the colored haplotypes
	my $pattern = $fields[2];
	for (my $strIdx = 0; $strIdx < $numStrains; $strIdx++) {
	    my $allele = substr($pattern,$strIdx, 1);
	    my $color = "white";
	    if ($allele ne "?") {
		$color = $haploColors[int($allele)];
	    }
	    print "<TD BGCOLOR=$color><FONT COLOR=$color SIZE=1>$allele</TD>\n";
	}
	print "</TD>";
	
	print "<TD>$fields[9]</TD>\n"; # chromosome name
	my $fpos = commify($fields[10]);
	my $lpos = commify($fields[11]);
	print "<TD>$fpos-$lpos</TD>\n"; # position of first-last SNPs in block

	if ($expressionFileFull) {
	    # show gene expression in C57BL6 tissues
	    for (my $i = 0; $i < $numtissues; $i++) {
		my $ch = substr($present, $i, 1);
		my $exprcolor;
		if (($ch eq 'A') || ($ch eq 'M')) {
		    $exprcolor = 'green';
		}
		elsif ($ch eq 'P') {
		    $exprcolor = 'blue';
		}
		elsif ($ch eq '-') {
		    $exprcolor = '#BEBEBE';
		}
		# FIXME:  Put in expresso link.
		print "<TD BGCOLOR=$exprcolor></TD>";
	    }
	}	
	print "</TR>\n";
    }
    close(RESULTS);
    print "</DIV>\n";
    print $tt->at_end;
    print p, end_html;
}

sub output_display
{
    print header;
    
    my $set_name = param("query_name");
    if ($set_name eq "") {
	$set_name = "Unnamed_data_set";
    }
    # untaint
    $set_name =~ /(.*)/;
    $set_name = $1;

    my $data_type = param("data_type");
    $data_type =~ /(Categorical|Quantitative)/;
    $data_type = $1;

    my $p_value = param("p_value");
    $p_value =~ /(.*)/;
    $p_value = $1;
    if ($p_value eq "") {
	$p_value = 0.05;
    }
    
    my $logreduce = param("logreduce");
    if (defined $logreduce) {
	$logreduce =~ /(.*)/;
	$logreduce = $1;
    }
    else {
	$logreduce = 'off';
    }

    my $codingOnly = param("codingOnly");
    if (defined $codingOnly) {
	$codingOnly =~ /(.*)/;
	$codingOnly = $1;
    }
    else {
	$codingOnly = 'off';
    }

    my $codingOnlyFlag = '';
    if ($codingOnly eq 'on') {
	$codingOnlyFlag = ' -f ';
    }

    my @regionChoices = param("regionChoices");
    
    my $nonOverlapping = param("nonOverlapping");
    if (defined $nonOverlapping) {
	$nonOverlapping =~ /(.*)/;
	$nonOverlapping = $1;
    }
    else {
	$nonOverlapping = 'off';
    }

 
    my $nonOverlappingFlag = '';
    if ($nonOverlapping eq 'on') {
	$nonOverlappingFlag = '-n';
    }
    
    my $goTerms = param("goTerms");

    # build table of phenotype values ($input_strains).
    my $input_strains = ();
    for my $strain (sort anumcmp keys %$strains) {
	$input_strains->{$strain} = param($strain);
    }

    # this is for the "force_eq" feature.
    my @equalClass = param("equalClass");


    $longnames = $longnames . "$nonOverlappingFlag";

# 	print " Set Name :$set_name\n";
# 	print "<br>";
# 	print " Data Type :$data_type\n";
# 	print "<br>";
# 	print " p-value :$p_value\n";
# 	print "<br>";

    # **** I think the real work starts here ****

    # build a unique file name for this combination of strains.
    # Concatenates all non-blank fields in order, then gets md5 digest
    for my $strain (sort anumcmp keys %$strains) {
	if ($input_strains->{$strain} ne "") {
	    $longnames = $longnames . "$strains->{$strain}";
	}
    }
    # print "longnames: $longnames<br>\n";

    my $digest = md5_hex($longnames);
    my $prefix = "$digest";

    my $rtmp = rand();
    my $formFile = md5_hex("$$" . "_" .  "$rtmp");

    my $eqFileFull = 0;
    if(@equalClass) {
	$eqFileFull = "$tmpdir/$prefix" . "_forceeq.txt";
	open(FORCE_EQ, ">$eqFileFull") || print "Failed to open $eqFileFull";    
	my %equalMap = map { $_ => 1 } @equalClass;
	for my $strain (sort anumcmp keys %$strains) {
## I think this is redundant
##	    $input_strains->{$strain} = param($strain);
	    next if ($input_strains->{$strain} eq "");
	    my $class = 0;
	    if (exists($equalMap{$strain})) {
		$class = 1;
	    }
	    print FORCE_EQ "$strain\t$class\n";
	    #write a zero
	}
	close FORCE_EQ;
    }

    (my $sec,my $usec) = gettimeofday();
    my $uniquePrefix = md5_hex("$sec$usec\n");#this is a prefix unique to this particular running of the hapmapper

#    my $strainsFileFull = "$dataDirPath/$prefix" . "_strains.txt";
    my $strainsFileFull = "$tmpdir/$prefix" . "_strains.txt";
    my $phenotypesFileFull = "$tmpdir/$uniquePrefix" . "_phenotypes.txt";
    my $haploblocksFileFull = "$tmpdir/$prefix" . "_haploblocks.txt";
    my $resultsFileFull = "$tmpdir/$uniquePrefix" . "_results.txt";
    my $resultsFileCpy; #we only use this in batch mode

    # results file needs to be saved for batch queries.
    if (defined $batchQuery) {
	$batchQuery =~ /.*[\/](.*)[.][^.]*/;
	my $batchPrefix = $1;
	$resultsFileCpy = "RESULTS/" . $batchPrefix . "_results.txt";
    }

    my $blockSNPsFile = "$tmpdir/$prefix" . "_SNPs.txt";
 
    open(STRAINS, ">$strainsFileFull") || print "<pre>Could not open for writing strains file $strainsFileFull: $!</pre>\n";
    open(PHENOTYPES, ">$phenotypesFileFull") || print "<pre>Could not open for writing phenotypes file $phenotypesFileFull: $!</pre>\n";
    # FIXME: check that you have a reasonable number of strains.
    for my $strain ( sort anumcmp keys %$strains) {
	$input_strains->{$strain} = param($strain);
	if ($input_strains->{$strain} ne "") {
	    # print the abbreviation to the STRAINS file
	    print STRAINS "$strain\t$strains->{$strain}\n";
	    my $strval = $input_strains->{$strain};
	    if ($logreduce eq 'on') {
		if ($strval <= 0) {
		    print "<center><b> Error: $strain has a value <= 0, which is illegal for log reduction.</b>\n";
		    print "<b> Please hit the \"back\" button and correct the values. </b></center>\n";
		    exit(1);
		}
		$strval = log($strval);
	    }
	    print PHENOTYPES "$strains->{$strain}\t$strval\n";
	}
    }
    close(PHENOTYPES);
    close(STRAINS);

    # Build haploblocks file if it does not yet exist for this combination of strains.
    # FIXME: user option to force rebuild.

#     FIXME:  Should check if SNPs file is there (also, whether minSNPs has changed!)
    if (-e $haploblocksFileFull) {
    }
    else {
	print "<center>There is no haplotype blocks file for this particular collection of strains.<br>\n";
	print "I will have to build one.  It should take less than ten minutes if there is no other<br>\n";
        print "heavy load on the system.  If it takes much longer than that, it is probably not going to finish.<br>\n";
        print "In that case, try hitting the \"back\" button and trying again, and, if that doesn't work,<br>\n";
        print "contact \"dill\@cs.stanford.edu\" for help.<p>\n";
	print "The next time you use the same set of strains, it should be significantly faster.<p>\n";
	print "</center>\n";

	opendir(CDIR, $SNPSdataDirPath) || print "<h3>Open of directory $SNPSdataDirPath failed: $!</h3>\n";
	my @chrfiles = grep(/.txt$/, readdir(CDIR));

	my @children = ();	# pids of child processes

	# FIXME: ? put a bound on number of simultaneous CPUs
	foreach my $file (@chrfiles) {
	    unless ($file eq "." || $file eq "..") {
		chomp($file);
		$file =~ /(\w+).txt/; # Get part of file name before ".txt".  It is the chromosome name.
		my $chrname = $1;
		# temporary file has process ID as suffix.  After successful build, this will be renamed.
		my $haploblocksFileTmp = $haploblocksFileFull . "_" . "$chrname" . "_" . "$$";
		my $blockSNPsFileTmp = $blockSNPsFile . "_" . "$chrname";

		my $pid = fork;
		if ($pid == 0) {
		    # I am the child
		    # close(STDOUT);
		    
		    my $eblockscmd = "$binDir/haplomap eblocks $nonOverlappingFlag -a $SNPSdataDirPath/$chrname.txt -g $geneCodingFile -s $strainsFileFull -o $haploblocksFileTmp -p $blockSNPsFileTmp";

		    #print "Launching $eblockscmd<br>\n";
		    {
			local %ENV = ();

			system($eblockscmd);
		    }
		    # FIXME: this needs to go in a log somewhere.
		    if ($? != 0) {
			print "<b>haplotype block construction failed with return code $?<br>Error $!<br>";
			print "Command: $eblockscmd<br>";
			print "This is a system error -- complain to dill\@cs.stanford.edu<b><br>\n";
			exit(1);
		    }
		    exit(0);	# if we don't exit here, child continues with loop!
		}
		else {
		    # I am the parent
		    # print "Spawned $pid<br>\n";
		    push @children, $pid;
		}
	    }
	}
	closedir(CDIR);

	# FIXME: Look into subtleties of waitpid options, return codes.
	foreach my $pid (@children) {
	    waitpid ($pid, 0);
	    # FIXME: check status code
	    # print "$pid finished<br>\n";
	}

	# Concatenate blocks for individual chromosomes into one file.
	my $catcmd = "cat $haploblocksFileFull" . "_" . '*' . "$$ > $haploblocksFileFull";
	my $rmcmd = "rm $haploblocksFileFull" . "_" . '*' . "$$";
	# print "Cat cmd: $catcmd<br>\n";
	{
	    local %ENV = ();
	    system($catcmd);
	}
	if ($?) {
	    print "Haplotype blocks file concatenation failed: $!<br>\n";
	    exit(1);
	}
	else {
	    # print "Haplotypes block file built.<br>\n";
	    {
		local %ENV = ();
		system($rmcmd);
	    }
	}

	# sort the haplotype blocks.
	my $sortblocksFile = "$haploblocksFileFull" . "_s_" . "$$";
	my $sortcmd = "perl sortblocks.pl $haploblocksFileFull $sortblocksFile";
	# print "Sorting cmd: $sortcmd\n";
	{
	    local %ENV = ();
	    system($sortcmd);
	}
	if ($?) {
	    print "Blocks file sorting failed: $!<br>\n";
	    exit(1);
	}
	else {
	    # print "Blocks file sorted.<br>\n";
	    if (!rename($sortblocksFile, $haploblocksFileFull)) {
		# FIXME:  Is there a way to get a specific error message from this?
		print "Rename of $sortblocksFile failed: $!\n";
		exit(1);
	    }
	    unlink("$sortblocksFile");
	}
    }

    my $dtarg = "";
    if ($data_type eq 'Categorical') {
	$dtarg = "-c ";
    }
    
    my $eqClassFlag = "";
    if($eqFileFull) {
	$eqClassFlag = "-q $eqFileFull ";
    }
    
    my $goFlag = "";
    my $goFileTmp = "$tmpdir/$prefix" . "_goterms.txt";
    if($goTerms) {
	open(GOFILE, ">$goFileTmp");
	print GOFILE "$goTerms\n";
	close(GOFILE);
	$goFlag = "-t $goFileFull -i $goFileTmp";
    }

    my $phmapcmd = ("$binDir/haplomap ghmap $dtarg -l $p_value -n '$set_name' " .
		    "-p $phenotypesFileFull -b $haploblocksFileFull " .
			((defined $popFileFull) ? " -r $popFileFull" : "") .
		    ((defined $expressionFileFull) ? " -e $expressionFileFull" : "") .
		    " $codingOnlyFlag -o $resultsFileFull $eqClassFlag $goFlag");

    # unlink old file, so we know if we failed
    if (-e $resultsFileFull) {
	unlink($resultsFileFull) || confess "Unlink of $resultsFileFull failed: $!\n";
    }

    # print "<h2>$phmapcmd</h2>\n";
    {
	local %ENV = ();
	system($phmapcmd);
    }
    # FIXME: It seems to be very difficult to capture an error message from system-invoked function.
    # see: http://perldoc.perl.org/perlfaq8.html    
     if (($? >> 8) != 0) {
 	print "<b>phenotype hapmapping failed with return code ", ($?>>8), "<br>Error $!<br>";
 	print "Command: $phmapcmd<br>\n";
 	print "This is a system error -- complain to dill\@cs.stanford.edu<b><br>\n";
 	exit(1);
     }

    results_html($data_type, $resultsFileFull, $prefix, $uniquePrefix,$p_value, @regionChoices);
    if (defined $batchQuery) {
	copy($resultsFileFull,$resultsFileCpy) || confess "could not copy results file to $resultsFileCpy";
    }
}


sub print_form
{
  
  print start_form;
  my $headergif = "";
  my @topheader = ();
  $headergif = td({-background =>"$images_url/body_bg.jpg", -align=>'CENTER', -height=>'120'},
		  big(font({-size=>6, -color=>'white'},"Haplo Map")));
  push( @topheader,$headergif);
  print center table ({-border => '0', -width => '100%'}, Tr(\@topheader));
  @topheader = ();
  my $rows = ();
  my $formrows = (); 
  
  print center "<Table width = 100%>";
  print "<td width= 20%>";
  print "</td>";
  print "<td bgcolor= '#FFFFFF', width= 60%>";  
  print "<br><br><br>";  
  $headergif = "";
  @topheader = ();
  $headergif = td({-background =>"$images_url/body_bg.jpg", -align=>'left'},
		  big(font({-size=>4, -color=>'white'},"Input data")));
  push( @topheader,$headergif);
  print center table ({-border => '0', -width => '100%'}, Tr(\@topheader)), p;
  @topheader = ();
  
#	$headergif = td({-bgcolor =>"#CFBD8B", -align=>'left'},big(font({-size=>4, -color=>'white'},"Enter values or Leave blank if no data")));
#	push( @topheader,$headergif);
#	print center table ({-border => '0', -width => '100%'}, Tr(\@topheader)), p;
#  @topheader = ();	
  
  my $count = 0;
  my $table_row;
  
  print "Enter request in form below or upload from local file: ";
  print "<input type=\"file\" name=\"uploaded_form\" />\n";
  print submit(-name => 'run', -value => 'Upload Form Fields'), p;

  push @$formrows, ( th ({-align=>'right'},"Set Name:") 
		     . td ({-align=>'left'},textfield (-name => 'query_name', -size => 20, -value => 'Unnamed_data_set'))
		     ); 
  
  push @$formrows, ( th ({-align=>'right'},"Data type:") 
		     . td ({-align=>'left'},popup_menu( -name => 'data_type',
							-value => ['Quantitative','Categorical'],
							-default => 'Quantitative')) ); 
  
  push @$formrows, ( th ({-align=>'right'},"p-value:") 
		     . td ({-align=>'left'},textfield (-name => 'p_value', -size => 10, -value => 0.05))
		     );      
  
  print table ({-border => '0'}, Tr($formrows)),p ;

  # print center table ({-border => '0'}, p;
  print checkbox(-name=>'logreduce', -label=>'Log reduce');
  print checkbox(-name=>'codingOnly', -label=>'Show genes with coding SNPs only'),p;

  print "Leave phenotype field blank if there is no data.\n";
  
  $rows = ();	   
  $table_row = "";
  for my $strain( sort anumcmp keys %$strains) {
      $table_row = $table_row . td ({-align=>'right'},"$strains->{$strain}")
	  . td (textfield (-name => "$strain", -size => 10)) . td ("     ");
      $count = $count + 1;
      if ($count == int($count/4)*4){
	  #$table_row = "(" . "$table_row" . ")";
	  push @$rows, ( $table_row );
	  
	  $table_row = "";
      }
  }
  if (($count % 4) > 0) {
      # there is a non-empty row to print.
      push @$rows, ( $table_row );
  }

  print center table ({-border => '0'}, Tr($rows)),p;
  

  if (defined $expressionFileFull) {
      print 'Show only genes expressed in these tissues:';
      print center checkbox_group(-name=>'regionChoices',
				  -values=>["aorta","ascending colon","back skin","bone marrow","cerebellum","cerebral cortex","descending colon",
			"diaphragm","endocrine pancreas","exocrine pancreas","gonadal fat pad","heart","heart left atrium",
			"heart left atrium and heart right atrium","heart left ventricle","heart right atrium","heart right ventricle",
			"hippocampus","interscapular brown adipose tissue","kidney","limb muscle","liver","lung","mammary gland",
			"mesenteric adipose tissue","spleen","striatum","subcutaneous adipose tissue","thymus",
			"tongue","trachea","urinary bladder"],
				  -columns=>5),p;
    #   print center checkbox_group(-name=>'regionChoices',
	# 			  -values=>["B Cell Basal", "Cerebellum", "Chondrocyte", 
	# 				    "Hippocampus", "Kidney", "Liver", "Lung", "Prefrontal Cortex",
	# 				    "Quadricep", "Spinal Cord", "Spleen", "Striatum", "T Cell Basal"],
	# 			  -columns=>5),p;
      print "Choose which strains must have the same haplotype:";
      my @values = sort anumcmp keys%{$strains};
      print center checkbox_group(-name=>'equalClass',
				  -values=>\@values,
				  -columns=>5),p;
      print "If you would like to limit to only certain GO terms, please enter those GO terms here, one per line.<br />";
      print center textarea(-name=>'goTerms',
			    -rows=>10,
			    -columns=>60);
      print "<br />";
  }
  print "<br />";

  print checkbox(-name=>'nonOverlapping', -label=>'Require non-overlapping blocks'),p;

  print center;
  print submit(-name => 'run', -value => 'Find');
  print reset(-name => 'Clear');
  print submit(-name => 'run', -value=>'Save Query');
  print p;
  print "</td>";
  print "<td width= 20%>";
  print "</td>";
  print "</Table>";

  # tell it what data set it is using
#  print "<h2> Putting SNPdata in hidden field: $SNPdata </h2>\n";
  print "<input type=hidden name=SNPdata value=$SNPdata />";
  print p;

#  For query upload, need the submit command ("Find") case above to 
#  see that it's a file, read file (param in "file handle context"), 
#  parse it, and then display a form with filled-in fields (if not too hard), which
#  can then be submitted as usual.
#  Also need to be able to download form data in same format.
#  print end_form();

  print end_html;
}

