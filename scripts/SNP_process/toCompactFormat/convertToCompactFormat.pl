#!/usr/local/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common -I/usr/lib/perl5/site_perl/5.8.8

use strict;
use FileHandle;
use commonFunctions;

my $inputFileHeader = "/home/mzheng/NGS/combined_results/05022012/chr";
my $inputAnnotationFile = "/home/mzheng/NGS/combined_results/05022012/SNP_codon_v65.txt";
my @interestedStrains = (0..21);
my @strainNames = ("129P2", "129S1", "129S5", "A_J", "AKR", "B_C", "C3H", "C57BL6NJ", "CBA", "DBA", "LPJ", "NOD", "NZO", "B10", "FVB", "LGJ", "MAMy", "MRL", "NZB", "NZW", "SMJ", "SJL");
my $outputFileHeader = "/home/mzheng/NGS/combined_results/05022012/compact/chr";

my @interestedChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
		     "13", "14", "15", "16", "17", "18", "19", "X");
#my @interestedChr = ("X");

my ($i, $j, $chr, $pos, $ref, $alt, $temp, $temp1);

die "error: size of strains not matched\n" if scalar(@interestedStrains) != scalar(@strainNames);

#### read the SNP annotations
my %SNPgene = ();
my %dist = ();
my %anno = ();
getSNPcodonInfo($inputAnnotationFile, \%SNPgene, \%dist, \%anno);
my $doingValidation = 1;
if($doingValidation == 1)
{
    foreach $chr (sort keys %SNPgene)
    {
	foreach $pos (sort {$a<=>$b;} keys %{$SNPgene{$chr}})
	{
	    my @parts1 = split /\|/, $SNPgene{$chr}->{$pos};
	    my @parts2 = split /\|/, $dist{$chr}->{$pos};
	    my @parts3 = split /\|/, $anno{$chr}->{$pos};
	    if(scalar(@parts1) != scalar(@parts2) || scalar(@parts1) != scalar(@parts3))
	    {
		print "$chr\t$pos\t$SNPgene{$chr}->{$pos}\t$dist{$chr}->{$pos}\t$anno{$chr}->{$pos}\n";
	    }
	}
	
    }
    exit(0);
}

### read SNPs (in NIEHS format)
my @inputFileList = ();
for($i = 0; $i <= $#interestedChr; $i++)
{
    push @inputFileList, $inputFileHeader.$interestedChr[$i].".txt";
}

my %SNPLocalID = ();
my %SNPSSID = ();
my %SNPAC = ();
my %SNPPos = ();
my %SNPRefAllele = ();
my %SNPAltAllele = ();
my %SNPAlleles = ();
commonFunctions::readSNPFileNIEHSFormat(\@inputFileList, \@interestedStrains, \@interestedChr, \%SNPLocalID, \%SNPSSID, \%SNPAC, \%SNPPos, \%SNPRefAllele, \%SNPAltAllele, \%SNPAlleles);

foreach $chr (@interestedChr)
{
    $temp1 = "chr$chr";
    my $output = new FileHandle;
    $output->open("> $outputFileHeader$chr.txt") or die "$outputFileHeader$chr.txt can't be written\n";
    print $output 'C57BL/6J';
    for($i = 0; $i < scalar(@strainNames); $i++)
    {
	print $output "\t$strainNames[$i]";
    }
    print $output "\n";
    for($i = 0; $i < scalar(@{$SNPPos{$chr}}); $i++)
    {
	$pos = $SNPPos{$chr}->[$i];
	die "error: annotation file not matched for chr$chr at $pos\n" if(!exists($SNPgene{$temp1}->{$pos}));
	$ref = $SNPRefAllele{$chr}->[$i];
	$alt = $SNPAltAllele{$chr}->[$i];
	print $output "$SNPAC{$chr}->[$i]\t$chr\t$pos\t$ref\/$alt\t0";
	for($j = 0; $j < scalar(@strainNames); $j++)
	{
	    $temp = $SNPAlleles{$chr}->[$i]->[$j];
	    if($temp eq $ref)
	    {
		print $output "0";
	    }
	    else
	    {
		if($temp eq $alt)
		{
		    print $output "1";
		}
		else
		{
		    print $output "?";
		}
	    }
	}	
	print $output "\t$SNPgene{$temp1}->{$pos}\t$dist{$temp1}->{$pos}\t$anno{$temp1}->{$pos}\n";
    }
    close($output);
}

print "succeeded\n";
exit(1);
sub getSNPcodonInfo
{
    my ($inputFile, $gene, $dist, $anno) = @_;
    die "insufficient num of parameters for getSNPcodonInfo\n" if !$anno;
    %{$gene} = ();
    %{$dist} = ();
    %{$anno} = ();
    my $input = new FileHandle;
    $input->open("< $inputFile") or die "$inputFile can't be read\n";
    my ($line, $lineCount) = ("", 0);
    my @parts;
    while(1)
    {
	$line = <$input>;
	last if !$line;
	chomp($line); chop $line if $line =~ /\r/;
	next if !$line;
	$lineCount++;
	print "current line: $lineCount\n" if($lineCount % 1000000 == 0);

	@parts = split /\t/, $line;
	die "failed: this line in $inputFile in bad format:\n$line\n" if(scalar(@parts) != 6);
	if(!exists($gene->{$parts[1]}))
	{
	    my %temp1 = ();
	    $gene->{$parts[1]} = \%temp1;
	    my %temp2 = ();
	    $dist->{$parts[1]} = \%temp2;
	    my %temp3 = ();
	    $anno->{$parts[1]} = \%temp3;
	}
	$gene->{$parts[1]}->{$parts[2]} = $parts[3];
	$dist->{$parts[1]}->{$parts[2]} = $parts[4];
	$anno->{$parts[1]}->{$parts[2]} = $parts[5];
    }
    close($input);
}
