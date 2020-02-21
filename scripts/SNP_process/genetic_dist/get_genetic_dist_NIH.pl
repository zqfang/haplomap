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

my $outputFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/genetic_dist/dist_NIH.txt";

## Roche SNP input;
my @row_array;
my $dbh = DBI->connect('DBI:Oracle:PELTZ', 'peltz_prod','test')
                or die "Couldn't connect to database: " . DBI->errstr;

my $sql = "select * from b36_MSNPDB_strain2";
my $sth = $dbh->prepare($sql);
$sth->execute;

my ($i, $j, $strain, $snp_id, $allele, $key, $index1, $index2);
my %allAllele = ();
my %strainIndex = ("129/Sv" => 0, "A/HeJ" => 1, "A/J"  => 2, "AKR/J" => 3, "B10.D2-H2/oSnJ" => 4, "BALB/cByJ" => 5,
		   "BALB/cJ" => 6, "BUB" => 7, "C3H/HeJ" => 8, "C57BL/6J" => 9, "DBA/2J" => 10, "FVB" => 11,
		   "LGJ" => 12, "LPJ" => 13, "MRL/MpJ" => 14, "NZB/BlnJ" => 15, "NZW/LaC" => 16, "SMJ" => 17);
my @allStrain = sort keys %strainIndex;
my $numberOfStrain = scalar(@allStrain);

my $count = 0;

while( @row_array = $sth->fetchrow_array)
{
    $strain = $row_array[1];
    if(!exists($strainIndex{$strain}))
    {
	next;
    }
    $allele = $row_array[2];
    if(uc($allele) ne "A" && uc($allele) ne "C" && uc($allele) ne "G" && uc($allele) ne "T")
    {
	next;
    }
    $snp_id = $row_array[3];
    if(!exists($allAllele{$snp_id}))
    {
	my @temp = ();
	$allAllele{$snp_id} = \@temp;
    }
    $allAllele{$snp_id}[$strainIndex{$strain}] = $allele;
    $count++;

    if($count % 100000 == 0)
    {
	print "$count: $strain\t$snp_id\t$allele\n";
    }
}
$sth->finish;
$dbh->disconnect;

my %isSNP = ();
my $countSNP = 0;
my $countAllele;
foreach $snp_id (sort keys %allAllele)
{
    my $isSNP = 0;
    my $ref = "";
    $countAllele = 0;
    for($i = 0; $i < $numberOfStrain; $i++)
    {
	if($allAllele{$snp_id}[$i])
	{
	    $countAllele++;
	    if(!$ref)
	    {
		$ref = $allAllele{$snp_id}[$i];
	    }
	    else
	    {
		if($ref ne $allAllele{$snp_id}[$i])
		{
		    $isSNP = 1;
		}
	    }
	}
    }
    if($countAllele < 9)
    {
	$isSNP = 0;
    }
    $isSNP{$snp_id} = $isSNP;
    $countSNP += $isSNP;
}

foreach $snp_id (keys %isSNP)
{
    if($isSNP{$snp_id})
    {
	for($i = 0; $i < $numberOfStrain; $i++)
	{
	    print "$allAllele{$snp_id}[$i],";
	}
	print "\n";
    }
}

###get genetic dist;
my @diff = ();
my @same = ();
my @missing = ();
for($i = 0; $i < $numberOfStrain; $i++)
{
    my @temp = ();
    $diff[$i] = \@temp;
    my @temp1 = ();
    $same[$i] = \@temp1;
    my @temp2 = ();
    $missing[$i] = \@temp2;
}

for($i = 1; $i < $numberOfStrain; $i++)
{
    $index1 = $strainIndex{$allStrain[$i]};
    for($j = 0; $j < $i; $j++)
    {
	$index2 = $strainIndex{$allStrain[$j]};
	my $countSame = 0;
	my $countDiff = 0;
	my $countMissing = 0;
	foreach $snp_id (sort keys %allAllele)
	{
	    if(!$isSNP{$snp_id})
	    {
		next;
	    }
	    if($allAllele{$snp_id}[$index1] && $allAllele{$snp_id}[$index2])
	    {
		if($allAllele{$snp_id}[$index1] ne $allAllele{$snp_id}[$index2])
		{
		    $countDiff++;
		}
		else
		{
		    $countSame++;
		}
	    }
	    else
	    {
		$countMissing++;
	    }
	}
	$diff[$i][$j] = $countDiff;
	$diff[$j][$i] = $countDiff;
	$same[$i][$j] = $countSame;
	$same[$j][$i] = $countSame;
	$missing[$i][$j] = $countMissing;
	$missing[$j][$i] = $countMissing;
    }
}

my $output = new FileHandle;
$output->open("> $outputFile") or die "$outputFile can't be opened.\n";
print $output "strain";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "\t$allStrain[$i]";
}
print $output "\n";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "$allStrain[$i]";
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
    print $output "\t$allStrain[$i]";
}
print $output "\n";
for($i = 0; $i < $numberOfStrain; $i++)
{
    print $output "$allStrain[$i]";
    for($j = 0; $j < $numberOfStrain; $j++)
    {
	if($i == $j)
	{
	    print $output "\t0";
	}
	else
	{
	    print $output "\t".(int($diff[$i][$j] / ($same[$i][$j] + $diff[$i][$j]) * 100) / 100);
	}
    }
    print $output "\n";
}

close($output);
print "succeeded\n";
