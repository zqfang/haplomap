#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common

use strict;
use FileHandle;
use Time::Local;

##### filtering parameters
my $qualCutoffForSamtools = 50;
my $cutoffForPLscoreForNonHomozygousAlt = 20;

my @strains = ("129P2", "129S1", "129S5", "AKR", "A_J", "B10", "BTBR", "BUB", "B_C", "C3H", "C57BL10J", "C57BL6NJ", "C57BRcd", "C57LJ", "C58", "CBA", "CEJ", "DBA", "DBA1J", "FVB", "ILNJ", "KK", "LGJ", "LPJ", "MAMy", "MRL", "NOD", "NON", "NUJ", "NZB", "NZO", "NZW", "PJ", "PLJ", "RFJ", "RHJ", "RIIIS", "SEA", "SJL", "SMJ", "ST", "SWR"); ### strains included in the VCF file in the exact order

my $inputFile = "/home/mzheng/NGS/combined_results/20160524/combined_gVCF/recal.SNP.vcf";
my $outputFolder = "/home/mzheng/NGS/combined_results/20160524/";

my @parts;
my @parts1;
my @parts2;
my @parts3;
my @parts4;
my @hasAlt;

my @alts;
my ($i, $j, $s, $D, $line, $GTind, $PLind, $ref, $alt, $lineCount, $numOfAlt, $index, $minScore, $thisChr, $thisPos, $numGoodAlt, $theGoodAlt);

my ($totalVariant, $totalInterested, $numNonPass, $numINDEL, $numLowQual, $numMultiAlt, $numGoodSNP) = (0, 0, 0, 0, 0, 0, 0);

my @chr = (10..19, 1..9, "X", "Y", "MT");
my %chr = ();
my %allAlleles = ();
my %allAlternatives = ();
my %allReferences = ();
for($i = 0; $i < scalar(@chr); $i++)
{
    $chr{"$chr[$i]"} = 1;
    my %temp = ();
    $allAlleles{$chr[$i]} = \%temp;
    my %temp1 = ();
    $allAlternatives{$chr[$i]} = \%temp1;
    my %temp2 = ();
    $allReferences{$chr[$i]} = \%temp2;
}

my %nucleotide = ();
$nucleotide{"A"} = 1; $nucleotide{"C"} = 1; $nucleotide{"G"} = 1; $nucleotide{"T"} = 1;

my $time1 = time();
my $inputVCF = new FileHandle;
$inputVCF->open("< $inputFile") or die "$inputFile can't be read\n";

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

    $totalVariant++;

    @parts = split /\t/, $line;
    die "format not recognized\n" if (scalar(@parts) != scalar(@strains) + 9);
    next if !exists($chr{$parts[0]}); ### not an interested Chr;
    $thisChr = $parts[0];
    $thisPos = $parts[1];

    $totalInterested++;
    if($parts[6] ne "PASS")
    {### does not pass the GATK filter, which is for SNP only!
	$numNonPass++;
	next;
    }

    ### first, get the ref/alternative alleles
    $ref = $parts[3];
    if(length($ref) > 1)
    {### an INDEL
	$numINDEL++;
	next;
    }

    $alt = $parts[4];
    @alts = split(/,/, $alt);
    $numOfAlt = scalar(@alts);
    
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

    @hasAlt = (0) x scalar(@alts);
    my @alleles = (-1) x scalar(@strains);
    for($s = 0; $s < scalar(@strains); $s++)
    {## identify the allele for each strain
	@parts2 = split(/\:/, $parts[$s + 9]);
	if(!$parts2[$GTind] || !$parts2[$PLind])
	{
##	    print "content for GT or PL not found properly:\n$line\n";
	    $alleles[$s] = "NN";
	    next;
	}
	@parts4 = split(/\//, $parts2[$GTind]);
	die "improper GT: $parts[$GTind]\n$line\n" if scalar(@parts4) != 2;
	if($parts4[0] != $parts4[1])
	{## not a good SNP
	    $alleles[$s] = "NN";
	    next;
	}
	## now this is supposed to be a homogyzous call. Check whether the call is of good quality
	$index = ($parts4[0] + 1) * ($parts4[0] + 2) / 2 - 1;
	@parts3 = split(/,/, $parts2[$PLind]);
	die "improper PL found: $parts2[$PLind]\n$line\n" if scalar(@parts3) != ($numOfAlt + 1) * ($numOfAlt + 2) / 2;
	$minScore = 10000000000;
	for($i = 0; $i < scalar(@parts3); $i++)
	{
	    next if $i == $index;
	    if($parts3[$i] - $parts3[$index] < $minScore)
	    {
		$minScore = $parts3[$i] - $parts3[$index];
	    }
	}
	if($minScore >= $cutoffForPLscoreForNonHomozygousAlt)
	{## this is a good call
	    if($parts4[0] > 0)
	    {
		$alleles[$s] = "$alts[$parts4[0] - 1]$alts[$parts4[0] - 1]";
	    }else
	    {
		$alleles[$s] = "$ref$ref";
	    }
	    if($parts4[0] > 0)
	    {
		$hasAlt[$parts4[0] - 1] = 1;
	    }
	}else
	{
	    $alleles[$s] = "NN";
	}
    }

    $numGoodAlt = 0;
    $theGoodAlt = -1;
    for($s = 0; $s < scalar(@hasAlt); $s++)
    {
	if($hasAlt[$s] == 1)
	{
	    $numGoodAlt++;
	    $theGoodAlt = $s;
	}
    }
    if($numGoodAlt == 1)
    {
	die "internal error: the good Alt not found properly" if $theGoodAlt == -1;
	if(!exists($nucleotide{$alts[$theGoodAlt]}))
	{## not SNP
	    $numINDEL++;
	}
	## this is a good SNP across the strains. save it for output
	$allAlleles{$thisChr}->{$thisPos} = \@alleles;
	$allAlternatives{$thisChr}->{$thisPos} = $alts[$theGoodAlt];
	$allReferences{$thisChr}->{$thisPos} = $ref;
	$numGoodSNP++;
    }else
    {
	if($numGoodAlt == 0)
	{### not a good SNP. no reliable homo alternative found.
	    $numLowQual++;
	}else
	{## multiple good alternatives found;
	    $numMultiAlt++;
	}
    }
}
close($inputVCF);

foreach $s (sort(keys %allAlleles))
{
    my $output = new FileHandle;
    $output->open("> $outputFolder"."chr$s.txt") or die ("$outputFolder"."chr$s.txt can't be written\n");
    ### print header
    print $output "LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES";
    for($j = 0; $j < scalar(@strains); $j++)
    {
	print $output "\t$strains[$j]";
    }
    print $output "\n";
    foreach $i (sort{$a<=>$b;} (keys%{$allAlleles{$s}}))
    {
	$D = "SNP_$s"."_$i";
	print $output "$D\t$D\t$s\t$D\t$i\t.\t$allReferences{$s}->{$i}\/$allAlternatives{$s}->{$i}";
	for($j = 0; $j < scalar(@strains); $j++)
	{
	    print $output "\t$allAlleles{$s}->{$i}->[$j]";
	}
	print $output "\n";
    }
    close($output);
}

print "total: $totalVariant, interested: $totalInterested, low-qual: $numNonPass, indels: $numINDEL, lowQual: $numLowQual, multiAlt: $numMultiAlt, good: $numGoodSNP\n";



