#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common

use strict;
use FileHandle;
use commonFunctions;
use DBI;

#############################################
# ENVIRONMENT VARIABLES#
#############################################
$ENV {'ORACLE_BASE'} = '/export/app/oracle';
$ENV{'ORACLE_SID'} = 'lseq';
$ENV{'ORACLE_TERM'} = 'vt100';
$ENV {'ORACLE_HOME'} = '/export/app/oracle/product/11.1.0';

#### input & output for NIEHS;
#my $outputFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/genetic_dist/dist_NIEHS.txt";
#my $inputHeader = "/home/mzheng/mzheng-data/SNP_data_perlegen/b04/b04_Chr";
#my $inputEnd = "_genotype.txt";
#my @allChr = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13",
#	      "14", "15", "16", "17", "18", "19", "X");
#my @allStrain = (7, 8, 9, 10, 11, 12, 13, 14, 17, 20, 21);
#my @strainName = ("C57", "DBA", "A/J", "B/B", "C3H", "AKR", "FVB", "129S1", "NOD", "BTBR", "NZW", "KK");


#### input & output for Sanger;
#my $outputFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/genetic_dist/dist_sanger.txt";
#my $inputHeader = "/home/mzheng/PeltzLabData/SNP_data/Sanger_01142010/ftp.sanger.ac.uk/pub/mouse_genomes/REL-0912/SNPs/processed/chr_";
#my $inputEnd = ".txt";
#my @allChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
#	      "14", "15", "16", "17", "18", "19", "X", "Y");
#my @allStrain = (7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19);
#my @strainName = ("C57/BL6", "129P2", "129S1", "129S5", "A_J", "AKR", "BALB", "C3H", "C57BL/6N", "CBA", "DBA", "LP_J", "NOD", "NZO");

#### input & output for combined;
#my $outputFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/genetic_dist/dist_combined.txt";
#my $outputFile1 = "/home/mzheng/mzheng-data/perl_prog/SNP_process/genetic_dist/kinship_matrix_combined.txt";
#my $inputHeader = "/home/mzheng/PeltzLabData/SNP_data/Sanger_01142010/ftp.sanger.ac.uk/pub/mouse_genomes/REL-0912/SNPs/processed/combined_Chr";
#my $inputEnd = ".txt";
#my @allChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
#	      "14", "15", "16", "17", "18", "19", "X");
#my @allStrain = (7..24);
#my @strainName = ("C57/BL6", "129S1/SvImJ", "A/J", "AKR/J", "C3H/HeJ", "DBA/2J", "NOD/ShiLtJ", "129P2", "129S5", "BALB/cJ", "C57BL/6NJ", "CBA/J", "LP/J", "NZO/HiLtJ", "BALB/cByJ", "FVB/NJ", "BTBR_T+_tf/J", "NZW/LacJ", "KK/HlJ");

my $outputFile = "/home/mzheng/NGS/combined_results/01242012/genetic_dist.txt";
my $outputFile1 = "/home/mzheng/NGS/combined_results/01242012/kinship_matrix.txt";
my $inputHeader = "/home/mzheng/NGS/combined_results/01242012/chr";
my $inputEnd = ".txt";
my @allChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
	      "14", "15", "16", "17", "18", "19", "X");
my @allStrain = (0..7, 9..13, 17..25);
my @strainName = ("C57/BL6", "129P2", "129S1", "129S5", "AKR", "A_J", "BALB", "C3H", "C57BL", "CBA", "DBA", "LP_J", "NOD", "NZO", "B10", "FVB", "LGJ", "MAMy", "MRL", "NZB", "NZW", "SMJ", "SJL");

my $numberOfChr = scalar(@allChr);
my $numberOfStrain = scalar(@strainName);

my($i, $j, $f, $line, $allele, $snp_id, $inputFile, $isSNP, $alleleCount, $isGoodSNP, $count, $countSNP, $ref, $alt);
my @parts = ();

$countSNP = 0;

my @indexRef;
my @indexAlt;
my @indexMissing;

###get genetic dist;
my @diff = ();
my @same = ();
my @missing = ();
for($i = 0; $i < $numberOfStrain; $i++)
{
    my @temp = (0) x scalar($numberOfStrain);
    $diff[$i] = \@temp;
    my @temp1 = (0) x scalar($numberOfStrain);
    $same[$i] = \@temp1;
    my @temp2 = (0) x scalar($numberOfStrain);
    $missing[$i] = \@temp2;
}

for($f = 0; $f < $numberOfChr; $f++)
{
    $inputFile = $inputHeader.$allChr[$f].$inputEnd;
    my $input = new FileHandle;
    $input->open("< $inputFile") or die "$inputFile can't be opened.\n";
    print "current: $inputFile\n";
    $count = 0;
    while(1)
    {
	$line = <$input>;
	last if (!$line);
	print "$count\n" if $count % 100000 == 0;
	$count++;
	next if($count == 1);
	chomp($line); chop($line) if $line =~ '\r';
	@parts = split /\t/, $line;
	$ref = substr($parts[6], 0, 1);
	$alt = "";
	@indexRef = (0);
	@indexAlt = ();
	@indexMissing = ();
	$isSNP = 0;
	$alleleCount = 0;
	$isGoodSNP = 1;
	for($i = 0; $i <= $#allStrain; $i++)
	{
	    $allele = substr($parts[$allStrain[$i]], 0, 1);
	    if($allele eq "A" || $allele eq "T" || $allele eq "C" || $allele eq "G")
	    {
		$alleleCount++;
		if($allele eq $ref)
		{
		    push @indexRef, ($i + 1);
		}
		else
		{
		    if(!$alt)
		    {
			$alt = $allele;
			push @indexAlt, ($i + 1);
			$isSNP = 1;
		    }
		    else
		    {
			if($allele ne $alt)
			{
			    print "SNP $parts[0] not bi-allelic and ignored\n";
			    $isGoodSNP = 0;
			}
			else
			{
			    push @indexAlt, ($i + 1);
			    $isSNP = 1;
			}
		    }
		}
	    }
	    else
	    {
		push @indexMissing, ($i + 1);
	    }
	}
	if($alleleCount < 7 || !$isGoodSNP)
	{
	    $isSNP = 0;
	}
	next if (!$isSNP);
	$countSNP++;
	die "interval error: ref or alt not found properly 1\n" if (scalar(@indexRef) == 0 || scalar(@indexAlt) == 0);
	die "interval error: ref or alt not found properly 2\n" 
	    if (scalar(@indexRef) + scalar(@indexAlt) + scalar(@indexMissing) != scalar(@strainName));
	for($i = 0; $i < scalar(@indexRef); $i++)
	{
	    for($j = 0; $j < scalar(@indexAlt); $j++)
	    {
		$diff[$indexRef[$i]][$indexAlt[$j]]++;
		$diff[$indexAlt[$j]][$indexRef[$i]]++;
	    }
	}
	for($i = 0; $i < scalar(@indexRef); $i++)
	{
	    $same[$indexRef[$i]][$indexRef[$i]]++;
	}
	for($i = 1; $i < scalar(@indexRef); $i++)
	{
	    for($j = 0; $j < $i; $j++)
	    {
		$same[$indexRef[$i]][$indexRef[$j]]++;
		$same[$indexRef[$j]][$indexRef[$i]]++;
	    }
	}
	for($i = 0; $i < scalar(@indexAlt); $i++)
	{
	    $same[$indexAlt[$i]][$indexAlt[$i]]++;
	}
	for($i = 1; $i < scalar(@indexAlt); $i++)
	{
	    for($j = 0; $j < $i; $j++)
	    {
		$same[$indexAlt[$i]][$indexAlt[$j]]++;
		$same[$indexAlt[$j]][$indexAlt[$i]]++;
	    }
	}
	for($j = 0; $j < scalar(@indexMissing); $j++)
	{
	    $missing[$indexMissing[$j]][$indexMissing[$j]]++;
	}
	for($j = 0; $j < scalar(@indexMissing); $j++)
	{
	    for($i = 0; $i < scalar(@indexRef); $i++)
	    {
		$missing[$indexRef[$i]][$indexMissing[$j]]++;
		$missing[$indexMissing[$j]][$indexRef[$i]]++;
	    }
	}	
	for($j = 0; $j < scalar(@indexMissing); $j++)
	{
	    for($i = 0; $i < scalar(@indexAlt); $i++)
	    {
		$missing[$indexAlt[$i]][$indexMissing[$j]]++;
		$missing[$indexMissing[$j]][$indexAlt[$i]]++;
	    }
	}
	for($i = 1; $i < scalar(@indexMissing); $i++)
	{
	    for($j = 0; $j < $i; $j++)
	    {
		$missing[$indexMissing[$i]][$indexMissing[$j]]++;
		$missing[$indexMissing[$j]][$indexMissing[$i]]++;
	    }
	}
    }
    
    close($input);
}

for($i = 0; $i < $numberOfStrain; $i++)
{
    for($j = 0; $j < $numberOfStrain; $j++)
    {
	die "error: sum of different situations not equal to total"
	    if($diff[$i][$j] + $same[$i][$j] +$missing[$i][$j] != $countSNP);
    }
}



my $output1 = new FileHandle;
$output1->open("> $outputFile1") or die "$outputFile1 can't be opened.\n";
print $output1 "strain";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output1 "\t$strainName[$i]";
}
print $output1 "\n";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output1 "$strainName[$i]";
    for($j = 0; $j < $numberOfStrain; $j++)
    {
	if($i == $j)
	{
	    print $output1 "\t".(($same[$i][$j] + 0.5 * $missing[$i][$j]) / $countSNP);
	}
	else
	{
	    print $output1 "\t".(($same[$i][$j] + 0.5 * $missing[$i][$j]) / $countSNP);
	}
    }
    print $output1 "\n";
}
close($output1);


my $output = new FileHandle;
$output->open("> $outputFile") or die "$outputFile can't be opened.\n";
print $output "strain";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "\t$strainName[$i]";
}
print $output "\n";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "$strainName[$i]";
    for($j = 0; $j < $numberOfStrain; $j++)
    {
	print $output "\t$diff[$i][$j]\/$same[$i][$j]\/$missing[$i][$j]";
    }
    print $output "\n";
}

print $output "total SNPs investigated: ".$countSNP."\n";


print $output "\n\nstrain";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "\t$strainName[$i]";
}
print $output "\n";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "$strainName[$i]";
    for($j = 0; $j < $numberOfStrain; $j++)
    {
	if($i == $j)
	{
	    print $output "\t0";
	}
	else
	{
	    print $output "\t".int($diff[$i][$j] / ($same[$i][$j] + $diff[$i][$j]) * 100) / 100;
	}
    }
    print $output "\n";
}

close($output);
print "succeeded\n";
