
# haplomap
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)


## Usage

### 1. convert vcf to niehs format for eblocks

SNPs
```shell
build/bin/haplomap niehs -o ${HOME}/data/SNPS/chr18.txt ${HOME}/data/VCFs/chr18.vcf

# support stdin, but much slower
zcat ${HOME}/data/VCFs/chr18.vcf.gz | bin/haplomap niehs -o ${HOME}/data/SNPS/chr18.txt
```
Structural variants
```shell
build/bin/haplomap niehs -o ${HOME}/data/SNPS/chr18.sv.txt --type sv input.sv.vcf
```


### 2. find haploblocks

prepared SNP annotation file 

1. generate strain level gene annotation database (only run once), see here: 
[scripts/gene_annotation](../scripts/gene_annotation/README.md)
    - this step generates 3 files for next step.
      - AA_by_strains_chr*.pkl 
      - mm10_kgXref.txt 
      - mm10_knownGene.txt

2. run `annotateSNPs.py` for each case (test_strains.txt) to get strain specific SNP annotation.
```shell
python scripts/annotateSNPs.py test_strains.txt chr18.txt \
                    AA_by_strains_chr18.pkl mm10_kgXref.txt mm10_knownGene.txt
                    genes_coding.txt genes_coding_transcript.txt
```

3. find haploblocks
```shell
build/bin/haplomap eblocks -a ${HOME}/data/SNPS/chr18.txt \
                     -g ${HOME}/data/gene_coding.txt \
                     -s ${HOME}/TMPDATA/test_strains.txt \
                     -o ${HOME}/TMPDATA/test.SNPs.hb.txt
```

prepared SV ananotation file
  - run `structral variant calling pipeline` and annotate with VEP
  - run `annotateSV.py`, the output file (*_eblocks.txt) is used for eblocks input.

### 3. statistical testing with trait data

Statistical testing  
  - SNP
    ```shell
    build/bin/haplomap ghmap -p ${HOME}/data/test_traits.txt \
                      -b ${HOME}/TMPDATA/test.SNPs.hb.txt \
                      -o ${HOME}/TMPDATA/test.final.output.txt
    ```
  - Structural variant  
    ```shell
    # not -k is needed for sv
    build/bin/haplomap ghmap -k -p ${HOME}/data/test_traits.txt \
                      -b ${HOME}/TMPDATA/test.sv.hb.txt \
                      -o ${HOME}/TMPDATA/test.sv.output.txt
    ```


**Note:** strain order in (-p) should keep the same to the (-b). That's, eblocks (-s)

## Input
1. eblocks:
    - Strain file (-s): 
      - two column txt file: |- abbrev -|- fullname -|
      - see `test.strain.txt` in the example folder
    - Allele file (-a): NIEHS compact format (use subcmd `niehs` to convert vcf to niehs)
    - Gene Annotation (-g): 
      - format: <SNP_{chr}_{postion}>  <gene_name>  < SNP_cateogry> 
      - see above to prepare this file

2. ghmap:
    - Trait file (-p):  
        - two column txt file: |- abbrev -|- value -|, see `test.trait.txt` in the example folder
        - same `abbrev value` pair could be set multiple times, but need to keep same order to strain pattern of (-b). Example:
        ```$xslt
           129S1	18.2
           129S1	19.1
           129S1	14.3
           129S1	17.2
           A_J	19.3
           A_J	18.2
           AKR	22.1
           AKR	20.0
           AKR	24.6
           AKR	21.4
        ```
    - haploblocks (-b): eblocks output file
    - genetic relation (-r): optional file, could obtain from plink pca

## Output

1. eblocks:

| NO | Field | Explanation |
|--- | ---- | ------------ |
|0 |chrom | chromosome idx      |
|1 |Block | Block id            |
|2 |Start | Block start idx     |
|3 |Size  | Block size          |
|4 |ChrBeg| Chromosome begin idx |
|5 |ChrEnd| Chromosome end idx  |
|6 |Pattern | Haplotype pattern |
|7 |GeneName| Associated Gene   |
|8 |Codon Map | Whether Change Codon |

2. ghmap:
  * gene-oriented results file

| NO |Field | Explanation |
|---| ---- | ------------ |
|0 |GeneName     | Associated Gene     |
|1 |CodonFlag    | condonChange ? 1:0  |
|2 |Pattern      | Haplotype pattern   |
|3 |FStat/Pvalue | isCategorical ? Fstat : Pvalue |
|4 |Effect       | Genetic Effect ( Omega^2 )   |
|5 |FDR          | Benjamini Hochberg. If categorical, skip |
|6 |popPvalue    | Pillai’s Trace Pvalue |
|7 |popFDR       | Pillai’s Trace FDR |
|8 |popYes       | Pillai’s Trace FDR Rejection | 
|9 |Chrom        | Chromosome idx      |
|10 |ChrBeg      | HaploBlock begin idx|
|11 |ChrEnd      | HaploBlock end idx  |
|12 |GeneExprMap | Gene expression Map |

  * block-oriented result file

BlockID | (IGNORED) | blockStart | blockSize | ChrIdx | ChrStart | ChrEnd | Pattern | Fstat/Pval | CodingMap ...



**CodonFlag**
* -1: Non-codon change
* 0: Synonymous (not important)
* 1: missense/nonsense...
* 2: Splicing site change

