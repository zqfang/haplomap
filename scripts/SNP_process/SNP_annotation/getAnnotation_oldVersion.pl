#!/usr/local/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common -I/usr/lib/perl5/site_perl/5.8.8

use strict;
use FileHandle;
use commonFunctions;
use NUCLEOTIDE;
use Scalar::Util qw(weaken isweak);
use Bio::SeqIO;

#{
#    my $seqio  = Bio::SeqIO->new(-file => "/home/mzheng/mzheng-data/ensembl_data/EMBL/mouse/B37/v58/Mus_musculus.0.dat", -format => "EMBL");
#    exit(0);
#}

my $version = "v65";

#my $inputFileHeader = "/home/mzheng/PeltzLabData/SNP_data/Sanger_07212011/combined/08072011/chr";
#my $inputFileHeader = "/home/mzheng/PeltzLabData/lab_users/kbusch/combined2/chr";

#my $inputFileHeader = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_03142012/mapping_results/formatted_EMMA_SNPs/betterEqual_chr";
#my @interestedStrains = (0);
#my $outputFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_03142012/mapping_results/formatted_EMMA_SNPs/"."betterEqual_SNPs_$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_03142012/mapping_results/formatted_EMMA_SNPs/codon_change_betterEqualSNPs_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_03142012/mapping_results/formatted_EMMA_SNPs/"."geneInfo_$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_03142012/mapping_results/formatted_EMMA_SNPs/"."badGene_$version.txt";

#my $inputFileHeader = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_12122012/formatted_EMMA_SNPs/better_chr";
#my @interestedStrains = (0);
#my $outputFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_12122012/formatted_EMMA_SNPs/"."better_SNPs_$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_12122012/formatted_EMMA_SNPs/codon_change_betterSNPs_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_12122012/formatted_EMMA_SNPs/"."geneInfo_$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/new_hapmapper_12122012/formatted_EMMA_SNPs/"."badGene_$version.txt";

#my $type = "betterEqual";
#my $inputFileHeader = "/home/mzheng/mzheng-data/projects/mouse_genetics/hapmapper/phenotype_data/result_EMMA/formatted_EMMA_SNPs/".$type."_chr";
#my @interestedStrains = (0);
#my $outputFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/hapmapper/phenotype_data/result_EMMA/formatted_EMMA_SNPs/".$type."_SNPs_$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/hapmapper/phenotype_data/result_EMMA/formatted_EMMA_SNPs/codon_change_".$type."SNPs_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/hapmapper/phenotype_data/result_EMMA/formatted_EMMA_SNPs/".$type."geneInfo_$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/mzheng-data/projects/mouse_genetics/hapmapper/phenotype_data/result_EMMA/formatted_EMMA_SNPs/".$type."badGene_$version.txt";

#my $inputFileHeader = "/home/mzheng/NGS/combined_results/05022012/chr";
#my @interestedStrains = (0..21);
#my $outputFile = "/home/mzheng/NGS/combined_results/05022012/SNP_codon_"."$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/NGS/combined_results/05022012/codon_change_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/NGS/combined_results/05022012/gene_info_"."$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/NGS/combined_results/05022012/bad_gene_info_"."$version.txt";

#my $inputFileHeader = "/home/mzheng/NGS/combined_results/20121212/chr";
#my @interestedStrains = (0..28);
#my $outputFile = "/home/mzheng/NGS/combined_results/20121212/SNP_codon_"."$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/NGS/combined_results/20121212/codon_change_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/NGS/combined_results/20121212/gene_info_"."$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/NGS/combined_results/20121212/bad_gene_info_"."$version.txt";

my $inputFileHeader = "/home/mzheng/NGS/combined_results/20131216/chr";
my @interestedStrains = (0..38);
my $outputFile = "/home/mzheng/NGS/combined_results/20131216/SNP_codon_"."$version.txt";
my $outputCodonChangeFile = "/home/mzheng/NGS/combined_results/20131216/codon_change_"."$version.txt";
my $outputGeneInfoFile = "/home/mzheng/NGS/combined_results/20131216/gene_info_"."$version.txt";
my $outputBadGeneInfoFile = "/home/mzheng/NGS/combined_results/20131216/bad_gene_info_"."$version.txt";

#my $inputFileHeader = "/home/mzheng/PeltzLabData/SNP_data/Sanger_07212011/combined/08112011/chr";
#my @interestedStrains = (0..17);
#my $outputFile = "/home/mzheng/PeltzLabData/SNP_data/Sanger_07212011/combined/08112011/SNP_codon_Sanger_inhouse_"."$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/PeltzLabData/SNP_data/Sanger_07212011/combined/08112011/codon_change_Sanger_inhouse_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/PeltzLabData/SNP_data/Sanger_07212011/combined/08112011/gene_info_Sanger_inhouse_"."$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/PeltzLabData/SNP_data/Sanger_07212011/combined/08112011/bad_gene_info_Sanger_inhouse_"."$version.txt";

#my $outputFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/SNP_codon_allSanger_"."$version.txt";
#my $outputCodonChangeFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/codon_change_allSanger_"."$version.txt";
#my $outputGeneInfoFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/gene_info_allSanger_"."$version.txt";
#my $outputBadGeneInfoFile = "/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/bad_gene_info_allSanger_"."$version.txt";

my $inputGeneFileHeader = "/home/mzheng/mzheng-data/ensembl_data/EMBL/mouse/B37/$version/Mus_musculus.";
my $exclusionListFile = "/home/mzheng/mzheng-data/ensembl_data/EMBL/mouse/B37/$version/gene_exclusion_list.txt";
my $numberOfGeneFiles = 30;

my @interestedChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
		     "13", "14", "15", "16", "17", "18", "19", "X");
#my @interestedChr = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "12", "14", "15", "16", "17", "18", "19", "X");

my %chrIsInterested = ();
foreach my $iChr (@interestedChr)
{
    $chrIsInterested{$iChr} = 1;
}

my %exclusionList = ();
getExclusionList($exclusionListFile, \%exclusionList);

my ($f, $i, $j, $k, $key, $line, $chr, $inputFileName, $count, $tempPos, $tempRef, $tempLocalID, $tempSSID, $tempAccNum, $isSNP, $isGood, $insertPos, $tempAlternative);
my ($indexSG, $indexEG, $SNP_CDS_Pos, $CDSseq, $CDS_AA);
my @indexSM;
my @indexEM;
my @indexSC;
my @indexEC;
my @startLoc; ### for mRNA
my @endLoc; ### for mRNA
my @sLoc; ### for CDS
my @eLoc; ### for CDS
my ($currentFile, $currentChr, $contigStart, $contigEnd, $currentGene, $ensGene, $geneDescription, $start, $end, $strand);
my ($isGoodGene, $currentTranscript, $geneIsGood);
my ($altDNA, $SNPRefAA, $SNPAltAA, $SNPRefAASeq, $SNPAltAASeq, $SNPAAPos);
my ($temp, $temp1, $temp2, $minDist, $geneSymbol, $maxFlanking);
my ($countGood, $countBad, $countNear, $countNotFound, $countIncon, $countImprove, $countPotential) = (0, 0, 0, 0, 0, 0, 0);
my @parts;
my @temp;
my @values;

### read SNPs (in NIEHS format)
my @inputFileList = ();
for($i = 0; $i <= $#interestedChr; $i++)
{
    push @inputFileList, $inputFileHeader.$interestedChr[$i].".txt";
}

#@inputFileList = ($inputFileHeader.$interestedChr[5].".txt");
#@inputFileList = ("/home/mzheng/mzheng-data/perl_prog/SNP_process/SNP_annotation/temp_partial_combined_Chr1.txt");
#@interestedChr = ("1");

for($i = 0; $i <= $#inputFileList; $i++)
{
#    print $inputFileList[$i]."\n";
}

my %SNPLocalID = ();
my %SNPSSID = ();
my %SNPAC = ();
my %SNPPos = ();
my %SNPRefAllele = ();
my %SNPAltAllele = ();
my %SNPAlleles = ();
commonFunctions::readSNPFileNIEHSFormat(\@inputFileList, \@interestedStrains, \@interestedChr, \%SNPLocalID, \%SNPSSID, \%SNPAC, \%SNPPos, \%SNPRefAllele, \%SNPAltAllele, \%SNPAlleles);

my %SNPCodon = ();
my %SNPGene = ();
my %SNPGeneDist = ();

my %SNPCodonChange = ();
#### initialize SNP codon to intergenic;
foreach $key (sort keys %SNPLocalID)
{
    my @temp1 = ("intergenic") x scalar(@{$SNPLocalID{$key}});
    $SNPCodon{$key} = \@temp1;
    my @temp2 = ("NULL") x scalar(@{$SNPLocalID{$key}});
    $SNPGene{$key} = \@temp2;
    my @temp3 = ("-1") x scalar(@{$SNPLocalID{$key}});
    $SNPGeneDist{$key} = \@temp3;
}

###all possible primary tags:
### >CDS<, >STS<, >exon<, >gene<, >mRNA<, >misc_RNA<, >misc_feature<, >source<

### all possible chr:
##>1<, >10<, >11<, >12<, >13<, >14<, >15<, >16<, >17<, >18<, >19<, 
##>2<, >3<, >4<, >5<, >6<, >7<, >8<, >9<, >MT<, >NT_??????<, >X<, >Y<

### gene information saved for SNPs that is outside the gene region.
my %geneName =();
my %geneEnsName = ();
my %geneStart = ();
my %geneEnd = ();
my %geneStrand = ();
my %geneFile = ();
my %geneDescription = ();

my @tempSNPCodon = (); ### codon for a gene
my @tempSNPGene = ();
my @tempSNPGeneDist = ();

my @tempSNPCodonCDS = (); ### codon for a transcript. A gene could have multiple transcript.

my %notGoodGenes = ();

#### first, get gene information
for($i = 0; $i < $numberOfGeneFiles; $i++)
#foreach $i (22)
#$i = 22;
{
    my ($seqio, $seqobj);
    $currentFile = $inputGeneFileHeader.($i * 100).".dat";
    print("current file: $currentFile\n");

    $seqio  = Bio::SeqIO->new(-file => $currentFile, -format => "EMBL");

    $count = 0;

    while($seqobj = $seqio->next_seq())
    {
	$count++;
	$currentGene = "NOGENESPECIFIED";
	$ensGene = "";
	$currentTranscript = "";
	$isGoodGene = 0;
        # gets sequence as a string from sequence object
	my $seqstr = $seqobj->seq(); # actual sequence as a string
	$temp = $seqobj->accession_number();
	$temp =~ s/\s+//g;
	@parts = split /\:/, $temp;
	($currentChr, $contigStart, $contigEnd) = ($parts[2], $parts[3], $parts[4]);
	print "count: $count\tchr$parts[2]\tstart: $parts[3]\tend: $parts[4]\tstrand: $parts[5]\n"
	    if ($count % 100 == 0);

	### only look at chr 1-19 and X.
	next if($currentChr =~ /NT_/);
	next if($currentChr eq "MT" || $currentChr eq "Y");

	### ignore if the chr is not interested
	next if (!exists($chrIsInterested{$currentChr}) || $chrIsInterested{$currentChr} == 0);

	if(!exists($geneName{$currentChr}))
	{
	    my %temp1 = ();
	    $geneName{$currentChr} = \%temp1;
	    my %temp2 = ();
	    $geneStart{$currentChr} = \%temp2;
	    my %temp3 = ();
	    $geneEnd{$currentChr} = \%temp3;
	    my %temp4 = ();
	    $geneEnsName{$currentChr} = \%temp4;
	    my %temp5 =();
	    $geneStrand{$currentChr} = \%temp5;
	    my %temp6 = ();
	    $geneDescription{$currentChr} = \%temp6;
	}

	for my $feat_object ($seqobj->get_SeqFeatures)
	{
#	    print "primary tag: ", $feat_object->primary_tag, "\n";
	    if(uc($feat_object->primary_tag) eq "GENE")
	    {
		if($currentGene ne "NOGENESPECIFIED" && $isGoodGene)
		{
		    ## clean up the annotation: remove dulicate entries.
		    for($j = 0; $j <= $#tempSNPCodon; $j++)
		    {
			if($tempSNPCodon[$j] =~ '!')
			{
			    @parts = split /\!/, $tempSNPCodon[$j];
			    my %temp = ();
			    for($k = 0; $k <= $#parts; $k++)
			    {
				$temp{$parts[$k]} = 1;
			    }
			    $tempSNPCodon[$j] = "";
			    foreach $k (sort keys %temp)
			    {
				if($tempSNPCodon[$j])
				{
				    $tempSNPCodon[$j] .= ('!'.$k);
				}
				else
				{
				    $tempSNPCodon[$j] = $k;
				}
			    }
			}
			else
			{
			    if($tempSNPCodon[$j] eq "unspecified")
			    {
				$tempSNPCodon[$j] = "misc_genic";
			    }
			}
		    }
		    for($j = 0; $j <= $#tempSNPCodon; $j++)
		    {
			if($SNPCodon{$currentChr}->[$j + $indexSG] eq "intergenic")
			{
			    $SNPCodon{$currentChr}->[$j + $indexSG] = $tempSNPCodon[$j];
			    $SNPGene{$currentChr}->[$j + $indexSG] = $tempSNPGene[$j];
			    $SNPGeneDist{$currentChr}->[$j + $indexSG] = $tempSNPGeneDist[$j];
			}
			else
			{
			    $SNPCodon{$currentChr}->[$j + $indexSG] .= ('|'.$tempSNPCodon[$j]);
			    $SNPGene{$currentChr}->[$j + $indexSG] .= ('|'.$tempSNPGene[$j]);
			    $SNPGeneDist{$currentChr}->[$j + $indexSG] .= ('|'.$tempSNPGeneDist[$j]);
			}
#			print $SNPLocalID{$currentChr}->[$indexSG + $j], "\tchr$currentChr\t", $SNPPos{$currentChr}->[$indexSG + $j], "\t$tempSNPCodon[$j]\n";
		    }
		}
		$start = $feat_object->location->start() + $contigStart - 1;
		$end = $feat_object->location->end() + $contigStart - 1;
		$strand = $feat_object->strand();
		@values = $feat_object->get_tag_values("locus_tag");
		die "gene in $currentFile at chr $currentChr, pos: $start has multiple symbols" if scalar(@values) > 1;

		$currentGene = $values[0];

		@values = $feat_object->get_tag_values("note");
		if(scalar(@values) == 0)
		{
		    $geneDescription = "NULL";
		}
		else
		{
		    if(scalar(@values) == 1)
		    {
			$geneDescription = $values[0];
		    }
		    else
		    {
			die "gene in $currentFile at chr $currentChr, pos: $start has multiple description";
		    }
		}
		@values = $feat_object->get_tag_values("gene");
		die "gene in $currentFile at chr $currentChr, pos: $start has multiple ensembl symbols" 
		    if scalar(@values) > 1;
		$ensGene = $values[0];

		$isGoodGene = 0; ### initialize gene as bad. it will be good when it gets the correct Amino Acid seq

		$geneIsGood = 1;
		###for genes in the exclusion list, we remove then;
		if(exists($exclusionList{$currentGene}))
		{
		    $geneIsGood = 0;
		    $notGoodGenes{$currentGene."_".$ensGene."_".$currentFile} = "excluded_manually";		    
		    next;
		}
		if($strand != 1 && $strand != -1)
		{
		    print ">$currentGene<\t>$start<\t>$end<\n";
#		    die "un-recognized strand: $strand, File $currentFile\n";
		    print "un-recognized strand: $strand, File $currentFile\n";
		    $geneIsGood = 0;
		    $notGoodGenes{$currentGene."_".$ensGene."_".$currentFile} = "unknown_strand";
		    next;
		}

		#### for genes looking like AC123456.7 or AC123456.78, we remove them
		if($currentGene =~ /[A-Z]{2}\d{6}.\d+/)
		{
		    $geneIsGood = 0;
		    $notGoodGenes{$currentGene."_".$ensGene."_".$currentFile} = "not-well-annotated";
		    next;
		}

		die "in $currentFile, gene $currentGene has multiple location, chr: $currentChr, pos: $start..$end\n" 
		    if $feat_object->location->isa('Bio::Location::SplitLocationI');
		if(exists($geneName{$currentChr}->{$start}))
		{
		    $geneName{$currentChr}->{$start} .= "\^".$currentGene;
		    $geneDescription{$currentChr}->{$start} .= "\^".$geneDescription;
		    $geneEnsName{$currentChr}->{$start} .= "\^".$ensGene;
		    $geneStart{$currentChr}->{$start} = $start;
		    $geneEnd{$currentChr}->{$start} .= "\^".$end;
		    $geneStrand{$currentChr}->{$start} .= "\^".$strand;
		    $geneFile{$currentChr}->{$start} .= "\^".$currentFile;
		}
		else
		{
		    $geneName{$currentChr}->{$start} = $currentGene;
		    $geneDescription{$currentChr}->{$start} = $geneDescription;
		    $geneEnsName{$currentChr}->{$start} = $ensGene;
		    $geneStart{$currentChr}->{$start} = $start;
		    $geneEnd{$currentChr}->{$start} = $end;
		    $geneStrand{$currentChr}->{$start} = $strand;
		    $geneFile{$currentChr}->{$start} = $currentFile;
		}
		
	#	print "$currentGene, $ensGene, $start, $end, $strand\n";

		if($start > $end)
		{
		    ($start, $end) = ($end, $start);
#		    print "gene $currentGene, file $currentFile, is in reverse order\n";
		}
		($indexSG, $indexEG) = getSNPIndex($start, $end, $SNPPos{$currentChr});
		if($indexEG > $indexSG)
		{
		    @tempSNPCodon = ("unspecified") x ($indexEG - $indexSG);
		    @tempSNPGene = ($currentGene) x ($indexEG - $indexSG);
		    @tempSNPGeneDist = (0) x ($indexEG - $indexSG);
		}
		else
		{
		    @tempSNPCodon = ();
		    @tempSNPGene = ();
		    @tempSNPGeneDist = ();
		}
	    }
	    if(uc($feat_object->primary_tag) eq "MRNA")
	    {
		if(!$geneIsGood)
		{
		    #### gene annotation is not good, so ignore this gene;
		    next;
		}
		#### match gene information and get transcript information.
		@temp = $feat_object->get_tag_values("gene");
		
		die "file $currentFile error: mRNA and Gene not match for $currentGene, $temp[0]\n" 
		    if $ensGene ne $temp[0];
		@temp = $feat_object->get_tag_values("note");
		@parts = split /=/, $temp[0];
		$currentTranscript = $parts[1];
		
		### get its location information
		@startLoc = ();
		@endLoc =();
		if ($feat_object->location->isa('Bio::Location::SplitLocationI'))
		{
		    my $index = 0;
		    for my $location ($feat_object->location->sub_Location())
		    {
			if($strand == 1)
			{
			    $startLoc[$index] = $location->start() + $contigStart - 1;
			    $endLoc[$index] = $location->end() + $contigStart - 1;
			}
			else
			{
			    unshift @startLoc, $location->start() + $contigStart - 1;
			    unshift @endLoc, $location->end() + $contigStart - 1;
			}
			$index++;
		    }
		}
		else
		{
		    $startLoc[0] = $feat_object->location->start() + $contigStart - 1;
		    $endLoc[0] = $feat_object->location->end() + $contigStart - 1;
		}

		if($indexEG > $indexSG)
		{
		    @tempSNPCodonCDS = ("unspecified") x ($indexEG - $indexSG);
		}
		else
		{
		    @tempSNPCodonCDS = ();
		}		

		@indexSM = (-1) x scalar(@startLoc);
		@indexEM = (-1) x scalar(@startLoc);
		for($j = 0; $j <= $#startLoc; $j++)
		{
		    ($indexSM[$j], $indexEM[$j]) = getSNPIndex($startLoc[$j], $endLoc[$j], $SNPPos{$currentChr});
		    die "interval error: mRNA is outside of gene in file $currentFile, gene $currentGene\n"
			if ($indexSM[$j] < $indexSG || $indexEM[$j] > $indexEG);
		    for(my $index = $indexSM[$j]; $index < $indexEM[$j]; $index++)
		    {
			$tempSNPCodonCDS[$index - $indexSG] = "EXON";
		    }
		}
		for($j = 0; $j < $#startLoc; $j++)
		{
		    for(my $index = $indexEM[$j]; $index < $indexSM[$j + 1]; $index++)
		    {
			die "internal error: index for SNPs in the intron not identified correctly"
			    if ($SNPPos{$currentChr}->[$index] <= $endLoc[$j] || $SNPPos{$currentChr}->[$index] >= $startLoc[$j + 1]);
			if($SNPPos{$currentChr}->[$index] - $endLoc[$j] <= 2 
			   || $startLoc[$j + 1] - $SNPPos{$currentChr}->[$index] <= 2)
			{
			    $tempSNPCodon[$index - $indexSG] = "SPLICE_SITE!INTRONIC";
			}
			else
			{
			    $tempSNPCodon[$index - $indexSG] = "INTRONIC";
			}
		    }
		}
		next;
	    }
	    if(uc($feat_object->primary_tag) eq "CDS")
	    {
		if(!$geneIsGood)
		{
		    #### gene annotation is not good, so ignore this gene;
		    next;
		}

		#### a transcript can only have one CDS
		#### match gene information and transcript information.
		@temp = $feat_object->get_tag_values("gene");
		die "file $currentFile error: CDS and Gene not match for $currentGene, $temp[0]\n" 
		    if $ensGene ne $temp[0];

		@temp = $feat_object->get_tag_values("note");
		@parts = split /=/, $temp[0];
		die "file $currentFile error: CDS and mRNA not match for $currentGene, $currentTranscript, $temp[0]\n" 
		    if $currentTranscript ne $parts[1];

		### get its location information and CDS DNA sequence
		@sLoc = ();
		@eLoc =();
		$CDSseq = "";
		if ($feat_object->location->isa('Bio::Location::SplitLocationI'))
		{
		    my $index = 0;
		    for my $location ($feat_object->location->sub_Location())
		    {
			if($strand == 1)
			{
			    $sLoc[$index] = $location->start() + $contigStart - 1;
			    $eLoc[$index] = $location->end() + $contigStart - 1;
			    $CDSseq .= substr($seqstr, $location->start() - 1, $location->end() - $location->start() + 1);
			}
			else
			{
			    unshift @sLoc, $location->start() + $contigStart - 1;
			    unshift @eLoc, $location->end() + $contigStart - 1;
			    $CDSseq .= commonFunctions::getReverseComplement(substr($seqstr, $location->start() - 1, $location->end() - $location->start() + 1));
			}
			$index++;
		    }
		}
		else
		{
		    $sLoc[0] = $feat_object->location->start() + $contigStart - 1;
		    $eLoc[0] = $feat_object->location->end() + $contigStart - 1;
		    if($strand == 1)
		    {
			$CDSseq = substr($seqstr, $feat_object->location->start() - 1, 
					 $feat_object->location->end() - $feat_object->location->start() + 1);
		    }
		    else
		    {
			$CDSseq = commonFunctions::getReverseComplement(substr($seqstr, $feat_object->location->start() - 1, $feat_object->location->end() - $feat_object->location->start() + 1));
		    }
		}
		$CDS_AA = NUCLEOTIDE::translate($CDSseq, 0);
		substr($CDS_AA, -1, 1) = "" if(substr($CDS_AA, -1, 1) eq "X" || substr($CDS_AA, -1, 1) eq '?');
		$CDS_AA =~ s/X/U/g;
		@temp = $feat_object->get_tag_values("translation");
		if(substr($temp[0], 0, 1) ne "M" || length($CDSseq) % 3 != 0)
		{
		    next;
		}
		if($CDS_AA ne $temp[0])
		{
		    my @temp1 = $feat_object->get_tag_values("protein_id");
#		    die "CDS AA seq not match for gene $currentGene in file $currentFile\n$CDS_AA(id: $temp1[0], length:",length($CDSseq),")\n$temp[0]\n" if (substr($CDS_AA, 0, 1) eq "M" && $CDS_AA ne $temp[0]);
		    if (substr($CDS_AA, 0, 1) eq "M" && $CDS_AA ne $temp[0])
		    {
			print "CDS AA seq not match for gene $currentGene in file $currentFile\n$CDS_AA(id: $temp1[0], length:",length($CDSseq),")\n$temp[0]\n";
			$notGoodGenes{$currentGene."_".$ensGene."_".$currentFile} = "CDS_not_match_".$CDS_AA."_$temp[0]";
			$geneIsGood = 0;
			next;
		    }
		    substr($temp[0], 0, 1) = "" if (substr($temp[0], 0, 1) eq "X" || substr($temp[0], 0, 1) eq "U");
		    my $t1 = NUCLEOTIDE::translate($CDSseq, 1);
		    $t1 =~ s/X/U/g;
		    substr($t1, -1, 1) = "" if (substr($t1, -1, 1) eq "U" || substr($t1, -1, 1) eq '?');
		    my $t2 = NUCLEOTIDE::translate($CDSseq, 2);
		    $t2 =~ s/X/U/g;
		    substr($t2, -1, 1) = "" if (substr($t2, -1, 1) eq "U" || substr($t2, -1, 1) eq '?');

		    if($temp[0] ne $CDS_AA && $temp[0] ne $t1 && $temp[0] ne $t2)
		    {
			print "$CDS_AA :+0\n";
			print "$t1 :+1\n";
			print "$t2 :+2\n";
			print "$temp[0]\n";
#			die "CDS AA seq not match for gene $currentGene in file $currentFile\n$CDS_AA(id: $temp1[0], length:",length($CDSseq),")\n$temp[0]\n";
			print "CDS AA seq not match for gene $currentGene in file $currentFile\n$CDS_AA(id: $temp1[0], length:",length($CDSseq),")\n$temp[0]\n";
			$notGoodGenes{$currentGene."_".$ensGene."_".$currentFile} = "CDS_not_match_".$CDS_AA."_$temp[0]";
			$geneIsGood = 0;
			next;		       
		    }
		}

###find SNPs that are in the CDS range and then find out whether it is a synonymous or non-synonymous coding		
		@indexSC = (-1) x scalar(@sLoc);
		@indexEC = (-1) x scalar(@sLoc);
		for($j = 0; $j <= $#sLoc; $j++)
		{
		    ($indexSC[$j], $indexEC[$j]) = getSNPIndex($sLoc[$j], $eLoc[$j], $SNPPos{$currentChr});
		    ### test CDS in mRNA?
		    die "interval error: CDS is outside of mRNA in file $currentFile, gene $currentGene\n"
			if ($indexSC[$j] < $indexSG || $indexEC[$j] > $indexEG);
		    for(my $index = $indexSC[$j]; $index < $indexEC[$j]; $index++)
		    {
			### check whether it is a synonymous or non-synonymous change
			$SNP_CDS_Pos = 0;
			$tempPos = $SNPPos{$currentChr}->[$index];	
			if($strand == 1)
			{
			    die "interval error, SNP $SNPLocalID{$currentChr}->[$index] is outside the range of CDS for gene $currentGene in file $currentFile\n"
				if($tempPos > $eLoc[$#eLoc]);
			    for($k = 0; $k <= $#sLoc; $k++)
			    {
				die "interval error, SNP $SNPLocalID{$currentChr}->[$index] is outside the range of mRNA for gene $currentGene in file $currentFile\n"
				    if($tempPos < $sLoc[$k]);
				if($tempPos <= $eLoc[$k])
				{
				    $SNP_CDS_Pos += ($tempPos - $sLoc[$k]);
				    last;
				}
				else
				{
				    $SNP_CDS_Pos += ($eLoc[$k] - $sLoc[$k] + 1);
				}
			    }
			}
			else
			{
			    die "interval error, SNP $SNPLocalID{$currentChr}->[$index] is outside the range of CDS for gene $currentGene in file $currentFile\n"
				if($tempPos < $sLoc[0]);
			    for($k = $#sLoc; $k >= 0; $k--)
			    {
				die "interval error, SNP $SNPLocalID{$currentChr}->[$index] is outside the range of mRNA for gene $currentGene in file $currentFile\n"
				    if($tempPos > $eLoc[$k]);
				if($tempPos >= $sLoc[$k])
				{
				    $SNP_CDS_Pos += ($eLoc[$k] - $tempPos);
				    last;
				}
				else
				{
				    $SNP_CDS_Pos += ($eLoc[$k] - $sLoc[$k] + 1);
				}
			    }			    
			}
			if($strand == 1)
			{
			    $altDNA = $SNPAltAllele{$currentChr}->[$index];
			}
			else
			{
			    $altDNA = commonFunctions::getReverseComplement($SNPAltAllele{$currentChr}->[$index]);
			}

			if($SNP_CDS_Pos % 3 == 0)
			{
			    die "internal error\n" if ($SNP_CDS_Pos + 2 > length($CDSseq));
			    $SNPRefAASeq = substr($CDSseq, $SNP_CDS_Pos, 3);
			    if($SNP_CDS_Pos / 3 < length($CDS_AA))
			    {
				## the stop codon
				$SNPRefAA = substr($CDS_AA, $SNP_CDS_Pos / 3, 1);
			    }else
			    {
				if($SNP_CDS_Pos / 3 == length($CDS_AA))
				{
				    $SNPRefAA = "X";
				}
				else
				{
				    die "SNP_CDS_Pos out of range: $currentGene\t$currentFile\n";
				}
			    }
			    $SNPAltAASeq = substr($CDSseq, $SNP_CDS_Pos, 3);
			    substr($SNPAltAASeq, 0, 1) = $altDNA;
			    $SNPAltAA = NUCLEOTIDE::translate($SNPAltAASeq, 0);
			    $SNPAAPos = $SNP_CDS_Pos / 3 + 1;
			}
			else
			{
			    if($SNP_CDS_Pos % 3 == 1)
			    {
				die "internal error\n" if ($SNP_CDS_Pos + 1 > length($CDSseq));
				$SNPRefAASeq = substr($CDSseq, $SNP_CDS_Pos - 1, 3);
				if(($SNP_CDS_Pos - 1) / 3 < length($CDS_AA))
				{
				    ## the stop codon
				    $SNPRefAA = substr($CDS_AA, ($SNP_CDS_Pos - 1) / 3, 1);
				}else
				{
				    if(($SNP_CDS_Pos - 1) / 3 == length($CDS_AA))
				    {
					$SNPRefAA = "X";
				    }
				    else
				    {
					die "SNP_CDS_Pos out of range: $currentGene\t$currentFile\n";
				    }
				}
				$SNPAltAASeq = substr($CDSseq, $SNP_CDS_Pos - 1, 3);
				substr($SNPAltAASeq, 1, 1) = $altDNA;
				$SNPAltAA = NUCLEOTIDE::translate($SNPAltAASeq, 0);
				$SNPAAPos = ($SNP_CDS_Pos - 1) / 3 + 1;
			    }
			    else ###$SNP_CDS_Pos % 3 == 2
			    {
				die "internal error\n" if ($SNP_CDS_Pos > length($CDSseq));
				$SNPRefAASeq = substr($CDSseq, $SNP_CDS_Pos - 2, 3);
				if(($SNP_CDS_Pos - 2) / 3 < length($CDS_AA))
				{
				    ## the stop codon
				    $SNPRefAA = substr($CDS_AA, ($SNP_CDS_Pos - 2) / 3, 1);
				}else
				{
				    if(($SNP_CDS_Pos - 2) / 3 == length($CDS_AA))
				    {
					$SNPRefAA = "X";
				    }
				    else
				    {
					die "SNP_CDS_Pos out of range: $currentGene\t$currentFile\n";
				    }
				}
				$SNPAltAASeq = substr($CDSseq, $SNP_CDS_Pos - 2, 3);
				substr($SNPAltAASeq, 2, 1) = $altDNA;
				$SNPAltAA = NUCLEOTIDE::translate($SNPAltAASeq, 0);
				$SNPAAPos = ($SNP_CDS_Pos - 2) / 3 + 1;
			    }
			}
			if($SNPRefAA eq $SNPAltAA)
			{
			    $tempSNPCodonCDS[$index - $indexSG] = "SYNONYMOUS_CODING";
			}
			else
			{
			    $tempSNPCodonCDS[$index - $indexSG] = "NON_SYNONYMOUS_CODING";
			}
			if(!exists($SNPCodonChange{$SNPLocalID{$currentChr}->[$index]}))
			{
			    $SNPCodonChange{$SNPLocalID{$currentChr}->[$index]} = "$currentGene:$SNPAAPos:$SNPRefAASeq/$SNPRefAA<->$SNPAltAASeq/$SNPAltAA";
			}
			else
			{
			    $SNPCodonChange{$SNPLocalID{$currentChr}->[$index]} .= "!$currentGene:$SNPAAPos:$SNPRefAASeq/$SNPRefAA<->$SNPAltAASeq/$SNPAltAA";
			}
		    }
		}
		### find 5prime and 3prime UTR here.
		my ($UTRs, $UTRe) = getSNPIndex($startLoc[0], $sLoc[0], $SNPPos{$currentChr});
		for(my $index = $UTRs; $index < $UTRe; $index++)
		{
		    if($SNPPos{$currentChr}->[$index] < $sLoc[0] && $tempSNPCodonCDS[$index - $indexSG] eq "EXON")
		    {
			if($strand == 1)
			{
			    $tempSNPCodonCDS[$index - $indexSG] = "5PRIME_UTR";
			}
			else
			{
			    $tempSNPCodonCDS[$index - $indexSG] = "3PRIME_UTR";
			}
		    }
		}
		($UTRs, $UTRe) = getSNPIndex($eLoc[$#eLoc], $endLoc[$#endLoc], $SNPPos{$currentChr});
		for(my $index = $UTRs; $index < $UTRe; $index++)
		{
		    if($SNPPos{$currentChr}->[$index] > $eLoc[$#eLoc] && $tempSNPCodonCDS[$index - $indexSG] eq "EXON")
		    {
			if($strand == 1)
			{
			    $tempSNPCodonCDS[$index - $indexSG] = "3PRIME_UTR";
			}
			else
			{
			    $tempSNPCodonCDS[$index - $indexSG] = "5PRIME_UTR";
			}			
		    }
		}
		#### verify that every SNP in the transcript range (not the gene range) is annotated properly.
		for($j = 0; $j <= $#tempSNPCodonCDS; $j++)
		{
		    if($tempSNPCodonCDS[$j] eq "EXON")
		    {
			for($k = 0; $k <= $#tempSNPCodonCDS; $k++)
			{
			    print $SNPLocalID{$currentChr}->[$indexSG + $k], "\t$currentChr\t", 
			    $SNPPos{$currentChr}->[$indexSG + $k] - $contigStart + 1, 
			    "\t$currentGene\t$tempSNPCodonCDS[$k]\n";
			}
			die "SNP for $currentGene, $currentTranscript in $currentFile not annotated properly\n";
		    }
		}
		$isGoodGene = 1;
		### now incorporate current mRNA into the corresponding gene
		for($j = 0; $j <= $#tempSNPCodonCDS; $j++)
		{
		    next if $tempSNPCodonCDS[$j] eq "unspecified";
		    if($tempSNPCodon[$j] eq "unspecified")
		    {
			$tempSNPCodon[$j] = $tempSNPCodonCDS[$j];
		    }
		    else
		    {
			$tempSNPCodon[$j] .= "\!$tempSNPCodonCDS[$j]";
		    }
		}		
		next;
	    }
	    if(uc($feat_object->primary_tag) eq "EXON")
	    {
		next;
	    }
	    if(uc($feat_object->primary_tag) eq "STS" || uc($feat_object->primary_tag) eq "MISC_RNA" 
	       || uc($feat_object->primary_tag) eq "MISC_FEATURE" || uc($feat_object->primary_tag) eq "SOURCE" )
	    {
		next;
	    }
	}
	### the last gene of the contig
	if($currentGene ne "NOGENESPECIFIED" && $isGoodGene)
	{
#	    print "gene symbol: $currentGene\n";
	    ## clean up the annotation: remove dulicate entries.
	    for($j = 0; $j <= $#tempSNPCodon; $j++)
	    {
		if($tempSNPCodon[$j] =~ '!')
		{
		    @parts = split /\!/, $tempSNPCodon[$j];
		    my %temp = ();
		    for($k = 0; $k <= $#parts; $k++)
		    {
			$temp{$parts[$k]} = 1;
		    }
		    $tempSNPCodon[$j] = "";
		    foreach $k (sort keys %temp)
		    {
			if($tempSNPCodon[$j])
			{
			    $tempSNPCodon[$j] .= ('!'.$k);
			}
			else
			{
			    $tempSNPCodon[$j] = $k;
			}
		    }
		}
		else
		{
		    if($tempSNPCodon[$j] eq "unspecified")
		    {
			$tempSNPCodon[$j] = "misc_genic";
		    }
		}
	    }
	    for($j = 0; $j <= $#tempSNPCodon; $j++)
	    {
		if($SNPCodon{$currentChr}->[$j + $indexSG] eq "intergenic")
		{
		    $SNPCodon{$currentChr}->[$j + $indexSG] = $tempSNPCodon[$j];
		    $SNPGene{$currentChr}->[$j + $indexSG] = $tempSNPGene[$j];
		    $SNPGeneDist{$currentChr}->[$j + $indexSG] = $tempSNPGeneDist[$j];
		}
		else
		{
		    $SNPCodon{$currentChr}->[$j + $indexSG] .= ('|'.$tempSNPCodon[$j]);
		    $SNPGene{$currentChr}->[$j + $indexSG] .= ('|'.$tempSNPGene[$j]);
		    $SNPGeneDist{$currentChr}->[$j + $indexSG] .= ('|'.$tempSNPGeneDist[$j]);
		}
#		print $SNPLocalID{$currentChr}->[$indexSG + $j], "\tchr$currentChr\t", $SNPPos{$currentChr}->[$indexSG + $j], "\t$tempSNPCodon[$j]\n";
	    }
	}
    }
}

#### finally for the intergenic SNPs, find the nearest gene and the distance
foreach $i (sort keys %SNPCodon)
{
    my @start = sort {$a <=> $b;} keys %{$geneStart{$i}};
    for($j = 0; $j < scalar(@{$SNPCodon{$i}}); $j++)
    {
	next if $SNPCodon{$i}->[$j] ne "intergenic";
	$tempPos = $SNPPos{$i}->[$j];
	my $ind = commonFunctions::getInsertPosition($tempPos, 1, \@start, scalar(@start));
	if($ind <= $#start)
	{
	    $minDist = $start[$ind] - $tempPos;
	    $temp = $geneName{$i}->{$start[$ind]};
	    if($temp =~ "\^")
	    {
		$temp =~ s/\^/\|/g;
	    }
	    $geneSymbol = $temp;
	    $maxFlanking = $minDist;
	}
	else
	{
	    die "internal error: empty chromosome found\n" if $ind == 0;
	    my $tempEnd = $geneEnd{$i}->{$start[$ind - 1]};
	    my $tempGene = $geneName{$i}->{$start[$ind - 1]};
	    if($tempEnd =~ "\^")
	    {
		my @p = split /\^/, $tempEnd;
		my @q = split /\^/, $tempGene;
		die "internal error: gene $tempGene and end $tempEnd not the same length\n" if scalar(@p) != scalar(@q);
		my $maxEnd = -1;
		my $maxGene = "";
		for(my $ll = 0; $ll < scalar(@p); $ll++)
		{
		    if($maxEnd < $p[$ll])
		    {
			$maxEnd = $p[$ll];
			$maxGene = $q[$ll];
		    }
		}
		die "interval error: maxEnd not found properly for $tempGene at end $tempEnd\n" 
		    if ($maxEnd == -1 || !$maxGene);
		$tempEnd = $maxEnd;
		$tempGene = $maxGene;
	    }
	    $minDist = $tempPos - $tempEnd;
	    if($minDist < 0)
	    {
		$minDist = 0;
	    }
	    $geneSymbol = $tempGene;
	    $maxFlanking = $minDist;
	}
	if($minDist > 0)
	{
	    for($k = $ind - 1; $k >= 0; $k--)
	    {
		last if ($tempPos - $start[$k] > 10000000 + $maxFlanking);
		my $tempEnd = $geneEnd{$i}->{$start[$k]};
		my $tempGene = $geneName{$i}->{$start[$k]};
		if($tempEnd =~ "\^")
		{
		    my @p = split /\^/, $tempEnd;
		    my @q = split /\^/, $tempGene;
		    die "internal error: gene $tempGene and end $tempEnd not the same length\n" 
			if scalar(@p) != scalar(@q);
		    my $maxEnd = -1;
		    my $maxGene = "";
		    for(my $ll = 0; $ll < scalar(@p); $ll++)
		    {
			if($maxEnd < $p[$ll])
			{
			    $maxEnd = $p[$ll];
			    $maxGene = $q[$ll];
			}
		    }
		    die "interval error: maxEnd not found properly for $tempGene at end $tempEnd\n" 
			if ($maxEnd == -1 || !$maxGene);
		    $tempEnd = $maxEnd;
		    $tempGene = $maxGene;
		}
		if($tempEnd >= $tempPos)
		{
		    $minDist = 0;
		    $geneSymbol = $tempGene;
		    last;
		}
		if($tempPos - $tempEnd < $minDist)
		{
		    $minDist = $tempPos - $tempEnd;
		    $geneSymbol = $tempGene;
		}
	    }
	}
	$SNPGene{$i}->[$j] = $geneSymbol;
	$SNPGeneDist{$i}->[$j] = $minDist;
	$SNPCodon{$i}->[$j] = "unannotated_genic" if $minDist == 0;
    }
}

my $output = new FileHandle;
$output->open(" >$outputFile") or die "$outputFile can't be opened\n";
foreach $i (sort keys %SNPCodon)
{
    for($j = 0; $j < scalar(@{$SNPCodon{$i}}); $j++)
    {
	print $output "$SNPLocalID{$i}->[$j]\tchr$i\t$SNPPos{$i}->[$j]\t$SNPGene{$i}->[$j]\t$SNPGeneDist{$i}->[$j]\t$SNPCodon{$i}->[$j]\n";
    }
}
close($output);

my $outputGeneInfo = new FileHandle;
$outputGeneInfo->open(" >$outputGeneInfoFile") or die "$outputGeneInfoFile can't be opened\n";
foreach $i (sort keys %geneName)
{
    foreach $j (sort {$a <=> $b;} keys  %{$geneStart{$i}})
    {
	print $outputGeneInfo "$geneName{$i}->{$j}\t$geneEnsName{$i}->{$j}\t$geneStart{$i}->{$j}\t$geneEnd{$i}->{$j}\t$geneStrand{$i}->{$j}\t$geneFile{$i}->{$j}\t$geneDescription{$i}->{$j}\n";
    }
}
close($outputGeneInfo);

my $outputBadGeneInfo = new FileHandle;
$outputBadGeneInfo->open(" >$outputBadGeneInfoFile") or die "$outputBadGeneInfoFile can't be opened\n";
foreach $i (sort keys %notGoodGenes)
{
    print $outputBadGeneInfo "$i\t$notGoodGenes{$i}\n";
}
close($outputBadGeneInfo);

my $outputCodonChange = new FileHandle;
$outputCodonChange->open("> $outputCodonChangeFile") or die "$outputCodonChangeFile can't be opened\n";
foreach $i (sort keys %SNPCodonChange)
{
    if($SNPCodonChange{$i} =~ '!')
    {
	@parts = split /\!/, $SNPCodonChange{$i};
	my %tt = ();
	for($j = 0; $j <= $#parts; $j++)
	{
	    $tt{$parts[$j]} = 1;
	}
	@temp = sort keys %tt;
	$temp = $temp[0];
	for($j = 1; $j <= $#temp; $j++)
	{
	    $temp .= '!'.$temp[$j];
	}
    }
    else
    {
	$temp = $SNPCodonChange{$i};
    }
    print $outputCodonChange "$i\t$temp\n";
}
close($outputCodonChange);


print "succeeded\n";
exit(1);


sub getSNPIndex
{
    ### for the array SNPPosArray, we have p0 < p1 < p2 < ....
    ### if p0 < p1 < $start < p2 <= $end < ...., then return (2, 3);
    ### if p0 < p1 = $start < p2 <= $end < ...., then return (1, 3);
    ### i.e. $start <= p_indexS < p2 <= $end < p_indexE
    my ($start, $end, $SNPPosArray) = @_;
    my ($indexS, $indexE);
    $indexS = commonFunctions::getInsertPosition($start, 1, $SNPPosArray, scalar(@{$SNPPosArray}));
    if($indexS > 0 && $start == $SNPPosArray->[$indexS - 1])
    {
	$indexS--;
    }
    $indexE = commonFunctions::getInsertPosition($end, 1, $SNPPosArray, scalar(@{$SNPPosArray}));
    die "internal error: end smaller than start\n" if $indexE < $indexS;
    return ($indexS, $indexE);
}

sub getExclusionList
{
    ### get the list of genes that should be excluded in the SNP annotation step
    ### this list is pre-determined.
    my ($inputFile, $list) = @_;
    die "insufficient num of input for getExclusionList\n" if(!$list);
    %{$list} = ();
    my $line;
    my $input = new FileHandle;
    $input->open("< $inputFile") or die "$inputFile can't be read\n";
    while(1)
    {
	$line = <$input>;
	last if !$line;
	chomp($line); chop($line) if $line =~ /\r/;
	next if !$line;
	$list->{$line} = 1;
    }
    close($input);
}
