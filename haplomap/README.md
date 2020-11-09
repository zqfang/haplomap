
# haplomap
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)


## Usage
Saome example to run `haplomap`

e.g.
```bash
# 1. convert vcf to niehs format for eblocks
bin/haplomap niehs -i ${HOME}/data/VCFs/chr18.vcf \
                   -o ${HOME}/data/SNPS/chr18.txt

# support stdin, but slower
zcat ${HOME}/data/VCFs/chr18.vcf.gz | bin/haplomap niehs -o ${HOME}/data/SNPS/chr18.txt

# 2. find haploblocks
bin/haplomap eblocks -a ${HOME}/data/SNPS/chr18.txt \
                     -g ${HOME}/data/gene_coding.txt \
                     -s ${HOME}/TMPDATA/test_strains.txt \
                     -o ${HOME}/TMPDATA/test.SNPs.hb.txt

# 3. statistical testing with trait data
bin/haplomap ghmap -p ${HOME}/data/test_traits.txt \
                   -b ${HOME}/TMPDATA/test.SNPs.hb.txt \
                   -o ${HOME}/TMPDATA/test.final.output.txt
```


## Input
1. eblocks:
    - Strain file (-s): column1 -> abbrev, column2 -> fullname
    - Allele file (-a): NIEHS compact format (use subcmd `niehs` to convert vcf to niehs)
    - Gene Annotation (-g): 
       - format: <SNP_{chr}_{postion}>  <gene_name>  < SNP_cateogry> 

2. ghmap:
    - Trait file (-p):  
        - tow column txt file: <abbrev> <value>
        - the <abbrev> <value> pair order should match to (-b).
        - same <abbrev> <value> pair could be set multiple times to input individual animal data. Example:
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

1. ebloks:

| NO | Field | Explaination |
|--- | ---- | ------------ |
|0 |chrom | chromosome idx      |
|1 |Block | Block id            |
|2 |Start | Block start idx     |
|3 |Size  | Block size          |
|4 |ChrBeg| Chrmosome begin idx |
|5 |ChrEnd| Chromosome end idx  |
|6 |Pattern | Haplotype pattern |
|7 |GeneName| Associated Gene   |
|8 |Codon Map | Whether Change Condon |

2. ghmap:
  * gene-oriented results file

| NO |Field | Explaination |
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

BlockID | (IGNORED) | blockStart | blockSize | ChrIdx | ChrStart | ChrEnd | Pattern | Fstat/Pval | CondingMap ...


## Changelog
v0.1
* A brand new cmdline interface `haplomap`, including sub-commands.
  - niehs: convert VCF to NIEHS compact format.
  - eblocks: find maximal haploblocks.
  - ghmap: haplotype statistical testing
  - pca: principal component analysis
  
* Refactor source code, Cross-platform support (Linux, MacOS)
  - replace all linux platform specific functions with STL
  - refactor to use modern C++11
  - fixed namespace contamination
  - re-design the whole source code, to make it more modularise and extensible 
  - add cmake support
  - add unitest support, much easier for debugging and development 

* ghmap:
  - cmdline usage improvement
  - add population structure testing (Pillai’s Trace)
  - add multiple hypothesis testing correction for annova and pillai's trace (Benjamini Hochberg procedure)
  - add support for raw animal data (individual data) input, which will give more statistical power..
  - fixed memory leak 
* eblocks:
  - cmdline usage improvement
  - fixed memory leak
  - make cmdline option (-p) become optional. 
* pca: 
  - add a new sub-command 
  - could be used for getting genetic relationship 
* niehs:
  - For historical reasons, eblocks use NIEHS compact format as input. 
  - To use haplomap more friendly, we now do this for you.  
 

v0.0
* The original version of haplomap