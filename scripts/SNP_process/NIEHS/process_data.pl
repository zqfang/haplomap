#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common



use strict;
use FileHandle;
use commonFunctions;

#### input path;
my $inputPath = "/home/mzheng/mzheng-data/SNP_data_perlegen/b04/";
#### output path for formated data
my $outputPath = "/home/mzheng/mzheng-data/SNP_data_perlegen/b04/processed/";

### input strain
### all strains are included.
my @interestedStrain = (0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 14);
my @interestedStrainName = ("DBA/2J", "A/J", "BALB/cByJ", "C3H/HeJ", "AKR/J", "FVB/NJ", "129S1/SvImJ", "NOD/LtJ",
			    "BTBR T+ tf/J", "NZW/LacJ", "KK/HlJ");

my @interestedChr = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15",
		     "16", "17", "18", "19", "X");

my $outputExtension = ".txt";

my ($f, $i, $line, $count, $tempAllele, $tempRef, $isSNP);
my @parts;

for($f = 0; $f <= $#interestedChr; $f++)
{
    ### read perlegen data
    my $inputFile1 = $inputPath."b04_Chr".$interestedChr[$f]."_genotype.txt";
    my $input1 = new FileHandle;
    $input1->open("< $inputFile1") or die "$inputFile1 can't be opened.\n";

    my $outputFile1 = $outputPath."b04_processed_Chr".$interestedChr[$f].$outputExtension;
    my $output1 = new FileHandle;
    $output1->open("> $outputFile1") or die "$outputFile1 can't be opened.\n";

    print $output1 "LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES";
    for($i = 0; $i <= $#interestedStrainName; $i++)
    {
	print $output1 "\t$interestedStrainName[$i]";
    }
    print $output1 "\n";

    $count = 0;

    while(1)
    {
	$line = <$input1>;
	last if(!$line);
	$count++;
	if($count == 1)
	{
	    next;
	}
	print "perlegen $interestedChr[$f] $count\n" if $count % 10000 == 0;
	chomp($line); chop($line) if $line =~ '\r';
	next if !$line;
	@parts = split /\t/, $line;
	$tempRef = substr($parts[6], 0, 1);
	$isSNP = 0;
	for($i = 0; $i <= $#interestedStrain; $i++)
	{
	    $tempAllele = substr($parts[7 + $interestedStrain[$i]], 0, 1);
	    if($tempAllele eq "A" || $tempAllele eq "C" || $tempAllele eq "G" || $tempAllele eq "T")
	    {
		if($tempAllele ne $tempRef)
		{
		    $isSNP = 1;
		    last;
		}
	    }
	}
	if($isSNP)
	{
	    print $output1 "$parts[0]\t$parts[1]\t$parts[2]\t$parts[3]\t$parts[4]\t$parts[5]\t$parts[6]";
	    for($i = 0; $i <= $#interestedStrain; $i++)
	    {
		print $output1 "\t$parts[7 + $interestedStrain[$i]]";
	    }
	    print $output1 "\n";
	}
    }
    close($input1);
    close($output1);
}

print "succeeded.\n";
