#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common

use strict;
use FileHandle;
use Time::Local;

##### filtering parameters
my $qualCutoff = 50;
my $cutoffForPLscoreForNonHomozygousAlt = 20;

## input SNPs for inbred strains
my $inputInbredHeader = "/home/mzheng/NGS/combined_results/20180101/noWildDerived/chr";
#my @chrs = ("MT");
my @chrs = (10..19, 1..9, "X", "Y", "MT");

## input wild-derived
my @wildStrains = ("CAST", "MOLF", "PWD", "PWK", "SPRET", "WSB");
#my @wildStrains = ("CAST");
my @wildDerivedVCF = (
    "/home/mzheng/NGS/CAST/BWA_results/SAM_results/sorted/output.raw.UG.vcf",
    "/home/mzheng/NGS/MOLF/BWA_results/SAM_results/sorted/output.raw.UG.vcf",
    "/home/mzheng/NGS/PWD/BWA_results/SAM_results/sorted/output.raw.UG.vcf",
    "/home/mzheng/NGS/PWK/BWA_results/SAM_results/sorted/output.raw.UG.vcf",
    "/home/mzheng/NGS/SPRET/BWA_results/SAM_results/sorted/output.raw.UG.vcf",
    "/home/mzheng/NGS/WSB/BWA_results/SAM_results/sorted/output.raw.UG.vcf");

my $outputFolder = "/home/mzheng/NGS/combined_results/20180101/";

my ($i, $line, $j, $chr, $pos, $lineCount, $header, $GTind, $PLind, $ref, $alt, $thisChr, $thisPos, $tempPos);
my @parts;
my @parts2;
my @parts3;


my %nucleotide = ();
$nucleotide{"A"} = 1; $nucleotide{"C"} = 1; $nucleotide{"G"} = 1; $nucleotide{"T"} = 1;

my %interestedChr = ();
my %inbredSNPs = ();
my %inbredRef = ();
my %inbredAlt = ();
my %wildDerivedSNPs = ();
for($i = 0; $i < scalar(@chrs); $i++)
{
    $interestedChr{"$chrs[$i]"} = 1;
    my %temp = ();
    $inbredSNPs{$chrs[$i]} = \%temp;
    my %temp1 = ();
    $inbredAlt{$chrs[$i]} = \%temp1;
    my %temp2 = ();
    $wildDerivedSNPs{$chrs[$i]} = \%temp2;
    my %temp4 = ();
    $inbredRef{$chrs[$i]} = \%temp4;
}
my $time1 = time();


for($i = 0; $i < scalar(@chrs); $i++)
{
    $chr = $chrs[$i];
    my $thisChr = new FileHandle;
    $thisChr->open("< $inputInbredHeader$chr.txt") or die ("< $inputInbredHeader$chr.txt can't be opened\n");

    $lineCount = 0;
    while(1)
    {
	$line = <$thisChr>;
	last if !$line;
	chomp($line); chop($line) if $line =~ /\r/;
	next if !$line;
	
	$lineCount++; print "line: $lineCount\n" if $lineCount % 100000 == 0;
	if($lineCount == 1)
	{
	    $header = $line;
	    next; ### ignore the title line
	}
	@parts = split /\t/, $line;
	$pos = $parts[4];
	$inbredSNPs{$chr}->{$pos} = $line;
	$inbredRef{$chr}->{$pos} = substr($parts[6], 0, 1).substr($parts[6], 0, 1);
	$inbredAlt{$chr}->{$pos} = substr($parts[6], 2, 1).substr($parts[6], 2, 1);
    }
    close($thisChr);
    foreach $pos (keys %{$inbredSNPs{$chr}})
    {
	my @temp3 = ($inbredRef{$chr}->{$pos}) x scalar(@wildDerivedVCF);
	$wildDerivedSNPs{$chr}->{$pos} = \@temp3;
    }
}

my @inconsistentCounts = (0) x scalar(@wildDerivedVCF);

for(my $wild = 0; $wild < scalar(@wildDerivedVCF); $wild++)
{
    my $inputVCF = new FileHandle;
    $inputVCF->open("< $wildDerivedVCF[$wild]") or die "$wildDerivedVCF[$wild] can't be read\n";
  
    $lineCount = 0;
    while(1)
    {
	$line = <$inputVCF>;
	last if !$line;
	chomp($line); chop($line) if $line =~ /\r/;
	next if !$line;
	
	next if(substr($line, 0, 1) eq "\#"); ### these are description lines;
	
	$lineCount++; print "line: $lineCount\n" if $lineCount % 100000 == 0;
	
#    last if $lineCount > 200000;
	@parts = split /\t/, $line;
	die "format not recognized\n" if (scalar(@parts) != 10); ### for single strain
	next if !exists($interestedChr{$parts[0]}); ### not an interested Chr;
	$thisChr = $parts[0];
	$thisPos = $parts[1];

	### first, get the ref/alternative alleles. For single strain, no need to consider multi-alternative case
	$ref = $parts[3];
	$alt = $parts[4];
	if(length($ref) > 1 || length($alt) > 1)
	{### an INDEL
	    for($j = 0; $j < length($ref); $j++)
	    {
		$tempPos = $thisPos + $j;
		if(exists($wildDerivedSNPs{$thisChr}->{$tempPos}))
		{
		    $wildDerivedSNPs{$thisChr}->{$tempPos}->[$wild] = "NN";
		}
	    }
	    next;
	}
	die "bad alt for this line\n" if !$nucleotide{$alt};

        ### INDEL cases dealt with already. So this variant is a SNP.
	if(!exists($wildDerivedSNPs{$thisChr}->{$thisPos}))
	{### this SNP is not in inbreds.
	    next;
	}

	if(($parts[6] ne "PASS" && $parts[6] ne ".") || $parts[5] < $qualCutoff)
	{### does not pass the GATK filter
	    $wildDerivedSNPs{$thisChr}->{$thisPos}->[$wild] = "NN";
	    next;
	}
	if($alt.$alt ne $inbredAlt{$thisChr}->{$thisPos})
	{## alt inconsistent;
	    $wildDerivedSNPs{$thisChr}->{$thisPos}->[$wild] = "NN";
	    next;
	}
	
###find the entries for GT & PL
	if($parts[8] eq 'GT:AD:DP:GQ:PL')
	{
	    $GTind = 0;
	    $PLind = 4;
	}else
	{
	    if($parts[8] eq 'GT:AD:DP:GQ:PGT:PID:PL')
	    {
		$GTind = 0;
		$PLind = 6;
	    }else
	    {
		@parts2 = split (/\:/, $parts[8]);
		$GTind = -1; $PLind = -1;
		for($i = 0; $i < scalar(@parts2); $i++)
		{
		    $GTind = $i if($parts2[$i] eq "GT");
		    $PLind = $i if($parts2[$i] eq "PL");
		}
		die "GT or PL not found for line:\n$line\n" if $GTind == -1 || $PLind == -1;
	    }
	}
	
	@parts2 = split(/\:/, $parts[9]);
	if(!$parts2[$GTind] || !$parts2[$PLind])
	{
	    print "content for GT or PL not found properly:\n$line\n";
#	    $alleles[$s] = "NN";
	    next;
	}
	@parts3 = split(/\//, $parts2[$GTind]);
	die "improper GT: $parts[$GTind]\n$line\n" if scalar(@parts3) != 2;
	if($parts3[0] != 1 || $parts3[1] != 1)
	{## not a good SNP
	    $wildDerivedSNPs{$thisChr}->{$thisPos}->[$wild] = "NN";
	    next;
	}

	## now this is supposed to be a homogyzous alternative call. Check whether the call is of good quality
	@parts3 = split(/,/, $parts2[$PLind]);
	die "improper PL found: $parts2[$PLind]\n$line\n" if scalar(@parts3) != 3;
	die "improper PL found: score for homo alt is not zero: $parts[$PLind]\n$line\n" if($parts3[2] != 0);
	if($parts3[0] >= $cutoffForPLscoreForNonHomozygousAlt && $parts3[1] >= $cutoffForPLscoreForNonHomozygousAlt)
	{### this is a good call
	    $wildDerivedSNPs{$thisChr}->{$thisPos}->[$wild] = "$alt$alt";
	}else
	{
	    $wildDerivedSNPs{$thisChr}->{$thisPos}->[$wild] = "NN";
	}
    }
    close($inputVCF); 
}

foreach $chr (sort {$a cmp $b;} keys %inbredSNPs)
{
    my $output = new FileHandle;
    $output->open("> $outputFolder"."chr$chr.txt") or die "$outputFolder"."chr$chr.txt can't be written\n";
    print $output $header;
    for($i = 0; $i < scalar(@wildStrains); $i++)
    {
	print $output "\t$wildStrains[$i]";
    }
    print $output "\n";
    foreach $pos (sort {$a <=> $b;} keys %{$inbredSNPs{$chr}})
    {
	print $output $inbredSNPs{$chr}->{$pos};
	for($i = 0; $i < scalar(@{$wildDerivedSNPs{$chr}->{$pos}}); $i++)
	{
	    print $output "\t$wildDerivedSNPs{$chr}->{$pos}->[$i]";
	}
	print $output "\n";
    }
    close($output);
}
 

exit(0);

