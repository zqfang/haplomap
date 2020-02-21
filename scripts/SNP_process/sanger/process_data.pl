#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common


#### this is to read and process the pileup data from Sanger's SNP data

use strict;
use FileHandle;
use commonFunctions;

my $includeCAST = 1;

#### input path;
my $inputPath = "/home/mzheng/PeltzLabData/SNP_data/Sanger_01142010/ftp.sanger.ac.uk/pub/mouse_genomes/REL-0912/SNPs/";
#### output path for formated data
my $outputFilePath = "/home/mzheng/PeltzLabData/SNP_data/Sanger_01142010/ftp.sanger.ac.uk/pub/mouse_genomes/REL-0912/SNPs/processed/";


### input strain
### the 4 wild-derived strains are excluded.
my @inputStrain = ("129P2", "129S1", "129S5", "A_J", "AKR", "BALB", "C3H", "C57BL", 
		   "CBA", "DBA", "LP_J", "NOD", "NZO");
my $outputExtension = ".txt";

if($includeCAST)
{
    @inputStrain = ("129P2", "129S1", "129S5", "A_J", "AKR", "BALB", "C3H", "C57BL", "CAST",
		    "CBA", "DBA", "LP_J", "NOD", "NZO");
    $outputExtension = "_withCAST.txt";    
}

### input strain
### all strains are included.
#my @inputStrain = ("129P2", "129S1", "129S5", "A_J", "AKR", "BALB", "C3H", "C57BL", 
#		   "CBA", "DBA", "LP_J", "NOD", "NZO", "CAST", "PWK", "SPRET", "WSB");
#my $outputExtension = "_full.txt";


my ($i, $j, $temp, $line, $key, $count, $type);
my @parts;


#### read input pileup file
my %SNP = ();
my %ref = ();
my %alternative = ();
my %problematicSNP = ();

for($i = 0; $i <= $#inputStrain; $i++)
{
    my $inputFile = $inputPath.$inputStrain[$i].".pileup";
    my $input = new FileHandle;
    $input->open("< $inputFile") or die "input file $inputFile can not be opened.\n";
    $count = 0;
    while(1)
    {
	$count++; print "$i\t$inputStrain[$i]\t$count\n" if ($count % 100000 == 0);
	$line = <$input>;
	last if !$line;
	chomp($line); chop $line if $line =~ '\r';
	next if !$line;
	@parts = split /\t/, $line;
	if($parts[4] < 10)
	{
	    next;
	}
	if($parts[2] eq "\*")
	{
	    next;
	}
	$key = $parts[0]."_".$parts[1];
	if(!exists($SNP{$key}))
	{
	    my @temp = ();
	    $SNP{$key} = \@temp;
	}
	if(!$ref{$key})
	{
	    $ref{$key} = $parts[2];
	}
	else
	{
	    if($ref{$key} ne $parts[2])
	    {
		$problematicSNP{$key} = 1;
	    }
	}
	if($parts[3] eq "A" || $parts[3] eq "C" || $parts[3] eq "G" || $parts[3] eq "T")
	{
	    $type = 1;
	}
	else
	{
	    if($parts[3] eq "W" || $parts[3] eq "S" || $parts[3] eq "M"
	       ||$parts[3] eq "K" || $parts[3] eq "R" || $parts[3] eq "Y")
	    {
		$type = 2;
	    }
	    else
	    {
		$type = 3;
	    }
	}
	if($type == 3)
	{
	    $problematicSNP{$key} = 4;
	    $SNP{$key}[$i] = "N";
	    next;
	}
	if(!$alternative{$key})
	{
	    $alternative{$key} = $parts[3];
	}
	else
	{
	    if($type == 2)
	    {
		if($alternative{$key} eq "A" || $alternative{$key} eq "C" || $alternative{$key} eq "G"
		   || $alternative{$key} eq "T")
		{
		}
		else
		{
		    if($alternative{$key} ne $parts[3])
		    {
			$problematicSNP{$key} = 2;
		    }
		}
	    }
	    else
	    {
		#type == 1 now
		if($alternative{$key} eq "A" || $alternative{$key} eq "C" || $alternative{$key} eq "G"
		   || $alternative{$key} eq "T")
		{
		    if($alternative{$key} ne $parts[3])
		    {
			$problematicSNP{$key} = 2;
		    }
		}
		else
		{
		    $alternative{$key} = $parts[3];
		}
	    }
	}
	$SNP{$key}[$i] = $parts[3];
    }
}

my %interestedChr = (
		     "1"=>1,
		     "2"=>1,
		     "3"=>1,
		     "4"=>1,
		     "5"=>1,
		     "6"=>1,
		     "7"=>1,
		     "8"=>1,
		     "9"=>1,
		     "10"=>1,
		     "11"=>1,
		     "12"=>1,
		     "13"=>1,
		     "14"=>1,
		     "15"=>1,
		     "16"=>1,
		     "17"=>1,
		     "18"=>1,
		     "19"=>1,
		     "X"=>1,
		     "Y"=>1
);

my @output  = ();
my @allChr = keys %interestedChr;
for($i = 0; $i < scalar(@allChr); $i++)
{
    $interestedChr{$allChr[$i]} = $i;
}
for($i = 0; $i < scalar(@allChr); $i++)
{
    $output[$i] = new FileHandle;
    $output[$i]->open("> $outputFilePath"."chr_".$allChr[$i].$outputExtension) or die "output $i can not be opened.\n";
    my $t = $output[$i];
    print $t "LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES";
    for($j = 0; $j <= $#inputStrain; $j++)
    {
	print $t "\t".$inputStrain[$j];
    }
    print $t "\n";
}

foreach $key (sort keys %SNP)
{
    if($problematicSNP{$key})
    {
	next;
    }
    @parts = split /_/, $key;
    if(!exists($interestedChr{$parts[0]}))
    {
	next;
    }
    if($alternative{$key} ne "A" && $alternative{$key} ne "C" && $alternative{$key} ne "G" && $alternative{$key} ne "T")
    {
	next;
    }
    my $t = $output[$interestedChr{$parts[0]}];
    print $t "SNP_"."$key\t"."SNP_"."$key\t$parts[0]\t"."SNP_"."$key\t"."$parts[1]\t.\t$ref{$key}"."/"."$alternative{$key}";
    for($i = 0; $i <= $#inputStrain; $i++)
    {
	if(!$SNP{$key}[$i])
	{
	    print $t "\t".$ref{$key}.$ref{$key};
	}
	else
	{
	    print $t "\t".$SNP{$key}[$i].$SNP{$key}[$i];
	}
    }
    print $t "\n";
}

for($i = 0; $i < scalar(@allChr); $i++)
{
    close($output[$i]);
}

my $outputProblematic = new FileHandle;
$outputProblematic->open("> $outputFilePath"."problematic.txt") or die "output $outputFilePath"."problematic.txt can not be opened.\n";
foreach $key (sort keys %problematicSNP)
{
    print $outputProblematic "$key\t$problematicSNP{$key}\n";
}
close($outputProblematic);


print "succeeded.\n";
