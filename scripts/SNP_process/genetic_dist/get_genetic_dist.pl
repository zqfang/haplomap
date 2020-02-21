#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common

use strict;
use FileHandle;
use commonFunctions;

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
my @strainName = ("129P2", "129S1", "129S5", "AKR", "A_J", "BALB", "C3H", "C57BL", "CBA", "DBA", "LP_J", "NOD", "NZO", "B10", "FVB", "LGJ", "MAMy", "MRL", "NZB", "NZW", "SMJ", "SJL", "C57/BL6");

my $numberOfChr = scalar(@allChr);
my $numberOfStrain = scalar(@strainName);

my($i, $j, $k, $f, $countSNP, $key, $temp, $temp1);

### input all data in NIEHS format
my @inputFileList = ();
for($i = 0; $i < scalar(@allChr); $i++)
{
    push @inputFileList, $inputHeader.$allChr[$i].$inputEnd;
}

my %SNPLocalID = ();
my %SNPSSID = ();
my %SNPAC = ();
my %SNPPos = ();
my %SNPRefAllele = ();
my %SNPAltAllele = ();
my %SNPAlleles = ();

commonFunctions::readSNPFileNIEHSFormat(\@inputFileList, \@allStrain, \@allChr, \%SNPLocalID, \%SNPSSID, \%SNPAC, \%SNPPos, \%SNPRefAllele, \%SNPAltAllele, \%SNPAlleles);

$countSNP = 0;

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

foreach $key (sort keys %SNPPos)
{
    for($i = 0; $i < scalar(@{$SNPAlleles{$key}}); $i++)
    {
	$countSNP++;
	for($j = 0; $j <= scalar(@{$SNPAlleles{$key}->[$i]}); $j++)
	{
	    if($j == scalar(@{$SNPAlleles{$key}->[$i]}))
	    {
		$temp = uc($SNPRefAllele{$key}->[$i]);
	    }
	    else
	    {
		$temp = uc($SNPAlleles{$key}->[$i]->[$j]);
	    }
	    if ($temp ne "A" && $temp ne "T" && $temp ne "C" && $temp ne "G")
	    {
		die "internal error: chr: $key, pos: $SNPPos{$key}->[$i], ref: $temp\n" if $j == scalar(@{$SNPAlleles{$key}->[$i]});
		$temp = "N";
		$missing[$j]->[$j]++;
	    }else
	    {
		$same[$j]->[$j]++;
	    }
	    
	    for($k = 0; $k < $j; $k++)
	    {
		$temp1 = uc($SNPAlleles{$key}->[$i]->[$k]);
		$temp1 = "N" if ($temp ne "A" && $temp ne "T" && $temp ne "C" && $temp ne "T");
		if($temp eq "N" || $temp1 eq "N")
		{
		    $missing[$j]->[$k]++;
		    $missing[$k]->[$j]++;
		}
		else
		{
		    if($temp eq $temp1)
		    {
			$same[$j]->[$k]++;
			$same[$k]->[$j]++;
		    }
		    else
		    {
			$diff[$j]->[$k]++;
			$diff[$k]->[$j]++;
		    }
		}
	    }
	}
    }
}


for($i = 0; $i < $numberOfStrain; $i++)
{
    for($j = 0; $j < $numberOfStrain; $j++)
    {
	die "error: sum of different situations not equal to total: $diff[$i][$j] \+ $same[$i][$j] \+ $missing[$i][$j] \!\= $countSNP\n"
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
	    print $output "\t".int($diff[$i][$j] / ($same[$i][$j] + $diff[$i][$j]) * 1000) / 1000;
	}
    }
    print $output "\n";
}

close($output);
print "succeeded\n";
