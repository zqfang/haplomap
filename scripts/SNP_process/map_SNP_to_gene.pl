#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common


#### this is to map SNP to its gene
#### default input format is perlegen's format.

use strict;
use FileHandle;
use commonFunctions;

my $flanking = 10000; #### 10K flanking range to call a SNP is in the neighborhood of a gene
my $maxLength = 2256181; ### the maximum length (in terms of genomic distance in bp) of a gene

#### gene info input
my $geneInfoFile = "/home/mzheng/mzheng-data/UCSCgenomeBrowser_data/mouse/B37/refFlat.txt";

#### SNP input
my $inputPath = "/home/mzheng/PeltzLabData/SNP_data/Sanger_01142010/ftp.sanger.ac.uk/pub/mouse_genomes/REL-0912/SNPs/processed/";

### SNP information output
my $outputFile = $inputPath."SNP_gene_info.txt";

my @interestedChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y");

my ($i, $j, $jj, $temp, $line, $key, $chr, $pos, $ID, $index, $inSymbol, $nearSymbol, $dist, $count);
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
    $parts[2] =~ s/chr//g;
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
    push @{$allChrSymbol{$parts[2]}}, $parts[0];
    push @{$allChrStart{$parts[2]}}, $parts[4];
    push @{$allChrEnd{$parts[2]}}, $parts[5];
    push @{$allChrStrain{$parts[2]}}, $parts[3];
}

#### sort the results by Start position
foreach $key (sort keys %allChrSymbol)
{
    for($i = 1; $i < scalar(@{$allChrStart{$key}}); $i++)
    {
	for($j = 0; $j < $i; $j++)
	{
	    if($allChrStart{$key}[$j] > $allChrStart{$key}[$i])
	    {
		my $tempSymbol = $allChrSymbol{$key}[$i];
		my $tempStart = $allChrStart{$key}[$i];
		my $tempEnd = $allChrEnd{$key}[$i];
		my $tempStrain = $allChrStrain{$key}[$i];
		for($jj = $i; $jj > $j; $jj--)
		{
		    $allChrSymbol{$key}[$jj] = $allChrSymbol{$key}[$jj - 1];
		    $allChrStart{$key}[$jj] = $allChrStart{$key}[$jj - 1];
		    $allChrEnd{$key}[$jj] = $allChrEnd{$key}[$jj - 1];
		    $allChrStrain{$key}[$jj] = $allChrStrain{$key}[$jj - 1];
		}
		$allChrSymbol{$key}[$j] = $tempSymbol;
		$allChrStart{$key}[$j] = $tempStart;
		$allChrEnd{$key}[$j] = $tempEnd;
		$allChrStrain{$key}[$j] = $tempStrain;
	    }
	}
    }
}

my $output = new FileHandle;
$output->open("> $outputFile") or die "output file $outputFile can't be opened\n";

for($j = 0; $j <= $#interestedChr; $j++)
#for($j = 20; $j < 21; $j++)
{
    my $inputFile = $inputPath."chr_".$interestedChr[$j].".txt";
    my $input = new FileHandle;
    $input->open("< $inputFile") or die "$inputFile can not be opened.\n";
    $count = 0;
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
	$ID = $parts[0];
	$chr = uc($parts[2]);
	$chr =~ s/CHR//g;
	$chr =~ s/_//g;
	$pos = $parts[4];
	if(!exists($allChrSymbol{$chr}))
	{
	    die "error: chr $chr not found\n";
	}
	$index = commonFunctions::getInsertPosition($pos, 1, $allChrStart{$chr}, scalar(@{$allChrStart{$chr}}));
	#### $allChrStart{$chr}[$index - 1] <= $pos < $allChrStart{$chr}[$index]
	$inSymbol = "NULL";
	$nearSymbol = "NULL";
	$dist = 10000000000;
	for($i = $index - 1; $i >= 0; $i--)
	{
	    if($pos - $allChrStart{$chr}[$i] > $flanking + $maxLength)
	    {
		last;
	    }
	    if($allChrEnd{$chr}[$i] >= $pos)
	    {
		$inSymbol = $allChrSymbol{$chr}[$i];
		$dist = 0;
		last;
	    }
	    if($pos - $allChrEnd{$chr}[$i] < $flanking)
	    {
		if($pos - $allChrEnd{$chr}[$i] < $dist)
		{
		    $dist = $pos - $allChrEnd{$chr}[$i];
		    $nearSymbol = $allChrSymbol{$chr}[$i];
		}
	    }
	}
	if($inSymbol eq "NULL" && $index < scalar(@{$allChrStart{$chr}}))
	{
	    if($allChrStart{$chr}[$index] - $pos < $flanking)
	    {
		if($allChrStart{$chr}[$index] - $pos < $dist)
		{
		    $dist = $allChrStart{$chr}[$index] - $pos;
		    $nearSymbol = $allChrSymbol{$chr}[$index];
		}
	    }
	}
	if($inSymbol ne "NULL")
	{
	    print $output "$ID\t$inSymbol\tIN\n";
	}
	else
	{
	    print $output "$ID\t$nearSymbol\tNEAR\n";
	}
    }
    close($input);
}

close($output);
print("succeeded.\n");
