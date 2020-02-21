#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common

#### check within a neighborhood of each gene, what is the SNP intensity.
#### default input format is perlegen's format.

use strict;
use FileHandle;
use commonFunctions;

my $flanking = 50000; #### 10K flanking range to call a SNP is in the neighborhood of a gene

#### gene info input
my $geneInfoFile = "/home/mzheng/mzheng-data/UCSCgenomeBrowser_data/mouse/B37/refFlat.txt";

#### SNP input
my $isAddingZero = 0;
my $inputPath = "/home/mzheng/PeltzLabData/SNP_data/Sanger_01142010/ftp.sanger.ac.uk/pub/mouse_genomes/REL-0912/SNPs/processed/";
my $inputHeader = "combined_Chr";

#my $isAddingZero = 1;
#my $inputPath = "/home/mzheng/mzheng-data/SNP_data_perlegen/b04/processed/";
#my $inputHeader = "b04_processed_Chr";

my @interestedChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X");

### SNP information output
my $outputFile = $inputPath."gene_SNP_coverage_info.txt";


my ($i, $j, $temp, $line, $key, $pos, $count, $s, $e, $indexS, $indexE, $isOld);
my @parts;
my @parts1;

##### read gene information
my $inputGeneInfo = new FileHandle;
$inputGeneInfo->open("< $geneInfoFile") or die "$geneInfoFile can't be opened\n";
my %allChrSymbol = ();
my %allChrStart = ();
my %allChrEnd = ();
my %allChrStrain = ();
while(1)
{
    $line = <$inputGeneInfo>;
    last if(!$line);
    chomp($line); chop($line) if $line =~ '\r';
    next if(!$line);
    @parts = split /\t/, $line;
    $parts[0] = uc($parts[0]);
    $parts[2] = uc($parts[2]);
    $parts[2] =~ s/CHR//g;
    if(!exists($allChrSymbol{$parts[2]}))
    {
	my @temp1 = ();
	$allChrSymbol{$parts[2]} = \@temp1;
	my @temp2 = ();
	$allChrStart{$parts[2]} = \@temp2;
	my @temp3 = ();
	$allChrEnd{$parts[2]} = \@temp3;
	my @temp4 = ();
	$allChrStrain{$parts[2]} = \@temp4;	
    }
    $isOld = 0;
    for($i = 0; $i < scalar(@{$allChrSymbol{$parts[2]}}); $i++)
    {
	if($allChrSymbol{$parts[2]}[$i] eq $parts[0])
	{
	    $isOld = 1;
	    last;
	}
    }
    if(!$isOld)
    {
	push @{$allChrSymbol{$parts[2]}}, $parts[0];
	push @{$allChrStart{$parts[2]}}, $parts[4];
	push @{$allChrEnd{$parts[2]}}, $parts[5];
	push @{$allChrStrain{$parts[2]}}, $parts[3];
    }
}

my $output = new FileHandle;
$output->open("> $outputFile") or die "output file $outputFile can't be opened\n";

print $output "chromosome\tsymbol\tgene_start\tgene_end\tflankingStart\tflankingEnd\tnumOfSNPinRegion\tSNPintensity(SNPperKB)\n";

for($j = 0; $j <= $#interestedChr; $j++)
#for($j = 20; $j < 21; $j++)
{
    my $inputFile;
    if($isAddingZero)
    {
	if($j < 9)
	{
	    $inputFile= $inputPath.$inputHeader."0".$interestedChr[$j].".txt";
	}
	else
	{
	    $inputFile= $inputPath.$inputHeader.$interestedChr[$j].".txt";
	}
    }
    else
    {
	$inputFile= $inputPath.$inputHeader.$interestedChr[$j].".txt";
    }
    my $input = new FileHandle;
    $input->open("< $inputFile") or die "$inputFile can not be opened.\n";
    $count = 0;
    my @allPos = ();
    while(1)
    {
	$line = <$input>;
	last if !$line;
	if($count == 0)
	{
	    $count++;
	    next;
	}
	$count++;
	chomp($line); chop($line) if $line =~ '\r';
	@parts = split /\t/, $line;
	$pos = $parts[4];
	push @allPos, $pos;
    }
    close($input);
    @allPos = sort {$a <=> $b} @allPos;
    
    if(!exists($allChrStart{$interestedChr[$j]}))
    {
	die "chromosome name inconsistent: $interestedChr[$j]\n";
    }
    for($i = 0; $i < scalar(@{$allChrStart{$interestedChr[$j]}}); $i++)
    {
	$s = $allChrStart{$interestedChr[$j]}[$i] - $flanking;
	$e = $allChrEnd{$interestedChr[$j]}[$i] + $flanking;	
	$indexS =commonFunctions::getInsertPosition($s, 1, \@allPos, scalar(@allPos));
	$indexE =commonFunctions::getInsertPosition($e, 1, \@allPos, scalar(@allPos));
	if($allPos[$indexS] == $s)
	{
	    $count = $indexE - $indexS + 1;
	}
	else
	{
	    $count = $indexE - $indexS;
	}
	print $output "chr_$interestedChr[$j]\t$allChrSymbol{$interestedChr[$j]}[$i]\t$allChrStart{$interestedChr[$j]}[$i]\t$allChrEnd{$interestedChr[$j]}[$i]\t$s\t$e\t$count\t".(($count / ($e - $s + 1)) * 1000)."\n";
    }
}

close($output);
print("succeeded.\n");
