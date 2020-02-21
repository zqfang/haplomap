#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common -I/home/mzheng/PeltzLabData/bioperl/bioperl-1.5.2_102/Bio

use strict;
use FileHandle;
use commonFunctions;

##### this is to compare the results with Sanger's annotation
#### first, read the sanger annotation.

my $inputSangerCodonFile = "/home/mzheng/mzheng-data/sanger_snps/SNP_annotations/REL_0912/snps-by-genes_13strains.tab";
#my $inputSangerCodonFile = "/home/mzheng/mzheng-data/sanger_snps/SNP_annotations/REL_0912/temp_snps_by_genes.txt";
my @interestedStrainsSanger = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);

my $inputCodonFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/SNP_codon_all_Sanger_v56.txt";

my $outputNotInSangerFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/not_in_Sanger.txt";
my $outputDiffGeneFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/SNP_with_diff_gene.txt";
my $outputDiffAnnoFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/SNP_with_diff_anno.txt";

my $outputDiffGenePairsFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/diff_gene_pairs.txt";

my ($i, $is, $j, $k, $tempChr, $key, $count, $line, $temp, $tempRef, $pos);
my @parts;

#### loading in house annotation
my $inputCodon = new FileHandle;
$inputCodon->open("< $inputCodonFile") or die "input file $inputCodonFile can't be opened\n";
my %codonGene1 = ();
my %codonType1 = ();
my %codonID = ();
my %codonDist = ();
$count = 0;
while(1)
{
    $line = <$inputCodon>;
    last if(!$line);
    chomp($line); chop($line) if $line =~ '\r';
    print "reading in house codon: $count\n" if $count % 100000 == 0;
    $count++;
    @parts = split /\t/, $line;
    $tempChr = uc($parts[1]);
    $tempChr =~ s/CHR//g;
    if(!exists($codonGene1{$tempChr}))
    {
	my %temp1 = ();
	$codonGene1{$tempChr} = \%temp1;
	my %temp2 = ();
	$codonType1{$tempChr} = \%temp2;
	my %temp3 = ();
	$codonID{$tempChr} = \%temp3;
	my %temp4 = ();
	$codonDist{$tempChr} = \%temp4;
    }
    $codonGene1{$tempChr}->{$parts[2]} = $parts[3];
    $codonType1{$tempChr}->{$parts[2]} = $parts[5];
    $codonID{$tempChr}->{$parts[2]} = $parts[0];
    $codonDist{$tempChr}->{$parts[2]} = $parts[4];
}
close($inputCodon);
#### loading Sanger's annotation

my $inputCodonS = new FileHandle;
$inputCodonS->open("< $inputSangerCodonFile") or die "input file $inputSangerCodonFile can't be opened.\n";
my %codonGene = ();
my %codonType = ();

$count = 0;

while(1)
{
    $line = <$inputCodonS>;
    last if(!$line);
    chomp($line); chop($line) if $line =~ '\r';
    print "reading Sanger codon: $count\n" if $count % 1000000 == 0;
    $count++;
    next if($count == 1);
    @parts = split /\t/, $line;
    $tempRef = $parts[3];

    my $isSNP = 0;
    my $isBad = 0;
    my $type = "";
    my $alternative = "";
    my $isBiallelic = 1;
    for($is = 0; $is <= $#interestedStrainsSanger; $is++)
    {
	$i = $interestedStrainsSanger[$is];
	if($parts[5 + $i * 2]) ###is different from the ref
	{
	    $temp = $parts[4 + $i * 2]; ###the alt allele
	    die "error for chr $parts[1], pos $parts[2]: allele >$temp<".
		"and annotation >".$parts[5 + $i * 2]."< inconsistent\n$line\n" if $temp eq "-";
	    $isSNP = 1;
	    if($temp eq "A" || $temp eq "C" || $temp eq "G" || $temp eq "T")
	    {
		if(!$alternative)
		{
		    $alternative = $temp;
		}
		else
		{
		    if($alternative ne $temp)
		    {
			$isBiallelic = 0;
			last;
		    }
		}
	    }
	    else
	    {
		#### alternvative looks like A/T
		my @parts3 = split /\//, $temp;
		if(scalar(@parts3) == 2 &&
		   (($parts3[0] eq "A" || $parts3[0] eq "C" || $parts3[0] eq "G" || $parts3[0] eq "T")
		    && ($parts3[1] eq "A" || $parts3[1] eq "C" || $parts3[1] eq "G" || $parts3[1] eq "T")))
		{
		    if($parts3[0] ne $tempRef)
		    {
			if(!$alternative)
			{
			    $alternative = $parts3[0];
			}
			else
			{
			    if($alternative ne $parts3[0])
			    {
				$isBiallelic = 0;
				last;
			    }
			}
		    }
		    if($parts3[1] ne $tempRef)
		    {
			if(!$alternative)
			{
			    $alternative = $parts3[1];
			}
			else
			{
			    if($alternative ne $parts3[1])
			    {
				$isBiallelic = 0;
				last;
			    }
			}
		    }
		}
		else
		{
		    ### un-recognized alternative
		    die("error1 for chr $parts[1], pos $parts[2]: irregular allele found, strain $i, allele: $temp\n");
		    $isBad = 1;
		    last;
		}
	    }
	    #### alternative allele checked now.

	    #### next, check the annotation.
	  #### in the case when alternative = A/T, one of them should be the reference, which doesn't have annotation
	    #### therefore, for annotation, don't need to consider the reference.
	    if(!$type)
	    {
		$type = $parts[5 + $i * 2];
	    }
	    else
	    {
		if($type ne $parts[5 + $i * 2])
		{
		    my @parts1 = sort(split /,/, $type);
		    my @parts2 = sort(split /,/, $parts[5 + $i * 2]);
		    my $isSame = 1;
		    if(scalar(@parts1) != scalar(@parts2))
		    {
			$isSame = 0;
		    }
		    else
		    {
			for($j = 0; $j <= $#parts1; $j++)
			{
			    if($parts1[$j] ne $parts2[$j])
			    {
				$isSame = 0;
				last;
			    }
			}
		    }
		    if(!$isSame)
		    {
			$isBad = 1;
			print "error for chr $parts[1], pos $parts[2]: annotation inconsistent for strains.\n";
		    }
		}
	    }
	}
	else
	{
	    die "error 2 for chr $parts[1], pos $parts[2]: allele and annotation inconsistent\n" 
		if $parts[4 + $i * 2] ne "-";
	}
    }
    
####save those good annotations
    if($isSNP & !$isBad & $isBiallelic)
    {
	if(!exists($codonGene{$parts[1]}))
	{
	    my %temp1 = ();
	    $codonGene{$parts[1]} = \%temp1;
	    my %temp2 = ();
	    $codonType{$parts[2]} = \%temp2;
	}
	if(exists($codonGene{$parts[1]}->{$parts[2]}))
	{
	    if($codonGene{$parts[1]}->{$parts[2]} eq $parts[0])
	    {
		die "chr $parts[1], pos $parts[2] already covered by the same gene and inconsistent.\n";
	    }
	    #### same SNP in a different genic region
	    $codonGene{$parts[1]}->{$parts[2]} .= ('|'.$parts[0]);
	    $codonType{$parts[1]}->{$parts[2]} .= ('|'.$type);
	}
	else
	{
	    $codonGene{$parts[1]}->{$parts[2]} = $parts[0];
	    $codonType{$parts[1]}->{$parts[2]} = $type;
	}
    }
}

close($inputCodonS);

my $outputNotInSanger = new FileHandle;
$outputNotInSanger->open("> $outputNotInSangerFile") or die "$outputNotInSangerFile can't be opened\n";

my $outputDiffGene = new FileHandle;
$outputDiffGene->open("> $outputDiffGeneFile") or die "$outputDiffGeneFile can't be opened\n";

my $outputDiffAnno = new FileHandle;
$outputDiffAnno->open("> $outputDiffAnnoFile") or die "$outputDiffAnnoFile can't be opened\n";

my $outputDiffGenePairs = new FileHandle;
$outputDiffGenePairs->open("> $outputDiffGenePairsFile") or die "$outputDiffGenePairsFile can't be opened\n";

my %diffGenePairs = ();

foreach $key (sort keys %codonGene1)
{
    if (!exists($codonGene{$key}))
    {
	print "chr $key not found in the sanger annotation\n";
	next;
    }
    print "processing chr$key...\n";
    $count = 0;
    foreach $pos (sort keys %{$codonGene1{$key}})
    {
	print "  processing SNP \# $count\n" if $count % 10000 == 0;
	$count++;
	if(!exists($codonGene{$key}->{$pos}))
	{
	    if($codonType1{$key}->{$pos} ne "intergenic" && $codonType1{$key}->{$pos} ne "misc_genic"
	       && $codonType1{$key}->{$pos} ne "unannotated_genic")
	    {
		print $outputNotInSanger "$codonID{$key}->{$pos}\tchr"."$key\t$pos\t$codonGene1{$key}->{$pos}\t$codonDist{$key}->{$pos}\t$codonType1{$key}->{$pos}\n";
	    }
	    next;
	}

	#### now the both annotations have this SNP

	### first,check the genes that is associated with the SNP
	my %SNPGeneTypeS = (); ### Sanger
	my %SNPGeneTypeI = (); ### in house;

	my @geneI = split /\|/, $codonGene1{$key}->{$pos};
	my @typeI = split /\|/, $codonType1{$key}->{$pos};
	die "annotation error (gene/type not matched) for $codonID{$key}->{$pos}\tchr"."$key\t$pos\t$codonGene1{$key}->{$pos}\t$codonDist{$key}->{$pos}\t$codonType1{$key}->{$pos}\n" if scalar(@geneI) != scalar(@typeI);
	for($j = 0; $j <= $#geneI; $j++)
	{
	    $SNPGeneTypeI{$geneI[$j]} = $typeI[$j];
	}



	my @geneS = split /\|/, $codonGene{$key}->{$pos};
	my @typeS = split /\|/, $codonType{$key}->{$pos};
	die "annotation error (gene/type not matched) for $codonID{$key}->{$pos}\tchr"."$key\t$pos\t$codonGene{$key}->{$pos}\t$codonType{$key}->{$pos}\n" if scalar(@geneS) != scalar(@typeS);
	for($j = 0; $j <= $#geneS; $j++)
	{
	    $SNPGeneTypeS{$geneS[$j]} = $typeS[$j];
	}
	
	my @keyS = sort(keys %SNPGeneTypeS);
	my @keyI = sort(keys %SNPGeneTypeI);
	if(scalar(@keyS) != scalar(@keyI))
	{
	    print $outputDiffGene "$codonID{$key}->{$pos}\tchr"."$key\t$pos\tinHouse:\t$codonGene1{$key}->{$pos}\t$codonDist{$key}->{$pos}\t$codonType1{$key}->{$pos}\tSanger:\t$codonGene{$key}->{$pos}\t$codonType{$key}->{$pos}\n";
	    $diffGenePairs{$codonGene1{$key}->{$pos}."_".$codonGene{$key}->{$pos}} = 1;
	    next;
	}
	my $genesAreSame = 1;
	for($j = 0; $j <= $#keyS; $j++)
	{
	    if($keyS[$j] ne $keyI[$j])
	    {
		$genesAreSame = 0;
		last;
	    }
	}
	if(!$genesAreSame)
	{
	    print $outputDiffGene "$codonID{$key}->{$pos}\tchr"."$key\t$pos\tinHouse:\t$codonGene1{$key}->{$pos}\t$codonDist{$key}->{$pos}\t$codonType1{$key}->{$pos}\tSanger:\t$codonGene{$key}->{$pos}\t$codonType{$key}->{$pos}\n";
	    $diffGenePairs{$codonGene1{$key}->{$pos}."_".$codonGene{$key}->{$pos}} = 1;
	    next;
	}

	### now, all genes are the same. Check each individual annotation term
	my $isSameAnnotation = 1;
	for($j = 0; $j <= $#keyS; $j++)
	{
	    my @typeI = sort(split /\!/, $SNPGeneTypeI{$keyI[$j]});
	    my @typeS = sort(split /\,/, $SNPGeneTypeS{$keyS[$j]});
	    if(scalar(@typeI) != scalar(@typeS))
	    {
		$isSameAnnotation = 0;
		last;
	    }
	    for($k = 0; $k <= $#typeI; $k++)
	    {
		if($typeI[$k] ne $typeS[$k])
		{
		    $isSameAnnotation = 0;
		    last;
		}
	    }
	    last if !$isSameAnnotation;
	}
	if(!$isSameAnnotation)
	{
	    print $outputDiffAnno "$codonID{$key}->{$pos}\tchr"."$key\t$pos\tinHouse:\t$codonGene1{$key}->{$pos}\t$codonDist{$key}->{$pos}\t$codonType1{$key}->{$pos}\tSanger:\t$codonGene{$key}->{$pos}\t$codonType{$key}->{$pos}\n";
	}
    }
}

foreach $key (sort keys %diffGenePairs)
{
    print $outputDiffGenePairs "$key\n";
}
close($outputNotInSangerFile);
close($outputDiffGene);
close($outputDiffAnno);
close($outputDiffGenePairs);
print "succeeded\n";
