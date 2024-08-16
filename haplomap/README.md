
# Haplomap
Haplotype-based computational genetic mapping  

## Installation

```shell
conda install -c bioconda haplomap
```

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)

![HBCGM](../docs/HBCGM.png)
## Usage
### 0. Variant calling
See [variant calling](../workflows/README.md) in the workflow folder using GATK, BCFtools, svtools.


Note: You need to use Ensemble-vep to annotate your VCF files.

the required flags: `--variant_class`, `--symbol`, `--individual_zyg all`, `--tab`

This code snap works for haplomap
```shell
## ensembl-vep version 101
bcftools view -f .,PASS ${vcf} | \
        vep --fasta ${reference} ${genome_build} \
        --format vcf --fork ${threads} --hgvs --force_overwrite \
        --uniprot --domains --symbol --regulatory --distance 1000 --biotype \
        --gene_phenotype MGI --check_existing  --pubmed --numbers \
        --offline --cache --variant_class \
        --gencode_basic --no_intergenic --individual all \
        -o ${output} --tab --compress_output gzip \

## for SV: https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#sv
## ensembl-vep version 111
bcftools view -f PASS ${vcf} | \
          vep -a GRCm38 --species mus_musculus --refseq --cache --cache_version 102 \
          --offline --compress_output gzip -o ${output} --tab \
          --fork 16 --offline --uniprot --cache --format vcf \
          --force_overwrite -overlaps --plugin TSSDistance --domains \
          --plugin phenotypes --symbol --canonical --variant_class \
          --nearest gene --regulatory --distance 5000 \
          --individual_zyg all --no_check_variants_order \
          --max_sv_size 100000 
```






### 1. Construct variant panel

**WARNING**: If you use GATK pipeline or Structrual Variants, you have to filter and subset SVs or Indels first. Haplomap `convert` work best with variants from BCFtools output. 

```shell
build/bin/haplomap convert --type snp \
                           -o ${HOME}/data/SNPS/chr18.txt ${HOME}/data/VCFs/chr18.vcf

# support stdin, but much slower
zcat ${HOME}/data/VCFs/chr18.vcf.gz | bin/haplomap convert -o ${HOME}/data/SNPS/chr18.txt
```




### 2. Construct annotation file from ensemble-vep results

NOTE:

1. Structural variants only support ensemble-vep inputs now !
2. Always check the variant name can be matched to `haplomap convert` output, since vep will change coordinate of variants in the output. Write your own annotation file if `haplomap annotate` did not work.

```shell
build/bin/haplomap annotate -o ${HOME}/data/SNPS/chr18.annotation.txt \
                            --type snp \
                            --samples test_strains.txt \ # only annotate the selected strains
                            input.vep.txt
```

**About file format**:

1. To annotate each variant, one gene follow by one annotation string (concat by ! with multiple consequence). 
2. For codon changes (missense, synonymous, stop), covert strings into a format like this: TCT/S<->CCT/P
3. To see all supported annotation, see [constants.cpp](src/constants.cpp)


e.g.
```
SNP_19_3358277  Cpt1a   splice_region_variant!AGG/R<->AGA/R
SNP_19_3365789  Cpt1a   splice_region_variant!intron_variant
SNP_19_3366324  Cpt1a   splice_region_variant!intron_variant
SNP_19_4000682  Nudt8   GAC/D<->GTC/V   Gm49405 GAC/D<->GTC/V!NMD_transcript_variant    Gm16312 downstream_gene_variant 4833408A19Rik   upstream_gene_variant
SNP_19_4000803  Nudt8   ACG/T<->ACT/T   Gm49405 ACG/T<->ACT/T!NMD_transcript_variant    Gm16312 non_coding_transcript_exon_variant      4833408A19Rik   upstream_gene_variant
SNP_19_4000809  Nudt8   CGT/R<->CGC/R   Gm49405 CGT/R<->CGC/R!NMD_transcript_variant    Gm16312 non_coding_transcript_exon_variant      4833408A19Rik   upstream_gene_variant 
SNP_19_4137099  Cabp4   splice_region_variant!intron_variant
SNP_19_4143385  Gpr152  CGA/R<->CGC/R!CGA/R<->CGT/R
```        


(Optional) If you'd like to use `ANNOVAR` (contributed by Boyoung Yoo @byoo1), see this

- 1. generate strain level gene annotation database (only run once), see here: 
[scripts/gene_annotation](../scripts/gene_annotation/README.md)
    - this step generates 3 files for next step.
      - AA_by_strains_chr*.pkl 
      - mm10_kgXref.txt 
      - mm10_knownGene.txt

- 2. run `annotateSNPs.py` for each case (test_strains.txt) to get strain specific SNP annotation.
```shell
python scripts/annotateSNPs.py test_strains.txt chr18.txt \
                    AA_by_strains_chr18.pkl mm10_kgXref.txt mm10_knownGene.txt
                    genes_coding.txt genes_coding_transcript.txt
```


### 3. Find haploblocks

```shell
build/bin/haplomap eblocks -a ${HOME}/data/SNPS/chr18.txt \
                     -g ${HOME}/data/chr18.annotation.txt \
                     -s ${HOME}/TMPDATA/test_strains.txt \
                     -o ${HOME}/TMPDATA/test.SNPs.hb.txt
```

### 4. Statistical testing with trait data

```shell
build/bin/haplomap ghmap -p ${HOME}/data/test_traits.txt \
                  -b ${HOME}/TMPDATA/test.SNPs.hb.txt \
                  -o ${HOME}/TMPDATA/test.final.output.txt
```
**NOTE 1**: 
By default. Output result are gene-summrized (see `haplomap ghmap --help`).   
Recommend adding `-a` flag, which will output gene-oriented format results.

**NOTE 2:** strain order in (-p) should keep the same to the eblocks (-s)


## Input
see `example` folder for test cases.

### convert
vcf file from variant calling tool.
vcf must be genotyped (has the GT field)

**Encoding**
1. SNV: A->A, T->T, G->G, C->C
2. INDEL/SV: 
   - reference: A
   - deletion: G
   - duplication or tanderm duplication: G
   - insertion: G
   - inversion: G
   - breakend: G


### annotate

For SVs, must contain SVTYPE, and SVLEN

### eblocks:
- Strain file (-s): 
  - Tree column txt file: "#Abbrev \t (Optional) \t Values "
  - see `test.strain.txt` in the example folder
- Allele file (-a): NIEHS compact format (use subcmd `convert` to convert vcf to niehs)
- Gene Annotation (-g): 
  - format: " <SNP_{chr}_{postion}>  <gene_name>  < consequence> "
  - see above to prepare this file

### ghmap:
- Trait file (-p):  
    - same as eblocks -s:  "#Abbrev \t (Optional) \t Values "
    - If multiple aninmal values for same strain, seperate them by comma. Example:
    ```$xslt
        129S1	18.2,19.1,14.3
        A_J	19.3,18.2
    ```
- haploblocks (-b): eblocks output file
- genetic distance matrix (-r): optional file, could obtain from plink.

**How to get the genetic distance matrix**

1. convert vcf to plink .tped, .tfam
```shell
haplomap convert  --plink \
                  -o ${HOME}/data/SNPS/chr1.snp.txt \
                  --type snp input.vcf
```

2. convert tped, tfam to .bed, .bim, .fam
```
plink --tfile chr1.snp --make-bed --out chr1
```
3. merge .bed files
```shell
plink --bfile chr1 --merge-list mergelist.txt --make-bed --out mouse_merged
```

note: the mergelist.txt in this format
  ```
  chr2.bed chr2.bim chr2.fam
  chr3.bed chr3.bim chr3.fam
  ...
  ```

4. Sample-distance and similarity matrices
```shell
plink --bfile mouse_merged # need .bim, .fam
      --make-rel square \ # Relationship/covariance
      --out mouse_grm \
      --threads 12
```
5. format to (ghmap -r) input
```shell
cut -f2 mouse_grm.rel.id | tr "\n" "\t" > mouse_grm.dist
echo "#$(cat mouse_grm.dist)" > mouse_grm.dist
cat mouse_grm.rel >> mouse_grm.dist
```


## Output

1. eblocks:

- Variant level file (-p): this file is used to check the individual variant in a haploblock

| NO | Field | Explanation |
|--- | ---- | ------------ |
|0 |Chrom | chromosome      |
|1 |Position | chromosome position (1-based) | 
|2 |SNP | SNP name |
|3 |Pattern | Binary pattern of alleles, strain order is same to the input strain file |
|4 |GeneName| Associated Gene   |
|5 |CodonFlag | Annotation |

- Aggregated haploblock file (-o): this file is used for ANNOVA

| NO | Field | Explanation |
|--- | ---- | ------------ |
|0 |Chrom | chromosome idx      |
|1 |BlockIdx | Block id            |
|2 |BlockStart | SNP vector start index; row index (0-based) of eblock -p output |
|3 |Size  | SNP number of the HaploBlock |
|4 |ChrBeg| Chromosome begin position (1-based) |
|5 |ChrEnd| Chromosome end position  |
|6 |Pattern | Haplotype pattern, strain order is same to the input strain file |
|7 |GeneName| Associated Gene   |
|8 |CodonFlag | Annotation |


2. ghmap:
  * gene-oriented results file

| NO |Field | Explanation |
|---| ---- | ------------ |
|0 |GeneName     | Associated Gene     |
|1 |CodonFlag    | [see here](src/constants.cpp)   |             
|2 |Haplotype    | Haplotype pattern, see header line for strain order   |
|3 |FStat/Pvalue | isCategorical ? Fstat : Pvalue |
|4 |EffectSize   | Genetic Effect ( Omega^2 )   |
|5 |FDR          | Benjamini Hochberg. If categorical, skip |
|6 |popPvalue    | Pillai’s Trace Pvalue |
|7 |popFDR       | Pillai’s Trace FDR |
|8 |Chr          | Chromosome      |
|9 |ChrStart     | Chromosome begin position (1-based) |
|10 |ChrEnd      | Chromosome end position   |
|11 | BlockIdx   | HaploBlock's index; the second column of eblock -o output |
|12 | BlockStart | SNPVector's index; row index (0-based) of eblock -p output |
|13 | BlockSize  | SNP number of HaploBlock |
|14 | Expression  | Gene expression Map, see header line  |

  * block-oriented result file

BlockIdx | BlockStart | blockSize | ChrIdx | ChrStart | ChrEnd | Pattern | Fstat/Pval | Effect | FDR | GeneExpMap | Gene | CodonFlag | ...



**CodonFlag**


i. SNPs
  * -1: Non-codon change
  * 0: Synonymous (not important)
  * 1: missense/nonsense
  * 2: Splicing site change
  * 3: Stop codon

ii. Indels and structral variants: 
  * 2: HIGH
  * 1: MODERATE
  * 0: LOW
  * -1: MODIFIER


See CodonFlag: [constants](src/constants.cpp) 

See explanation [here](https://ensembl.org/info/genome/variation/prediction/predicted_data.html) 

