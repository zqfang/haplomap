
# haplomap
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)


## Usage
e.g.
```bash
# find haplotypes
bin/haplomap eblocks -a ${HOME}/data/SNPS/chr18.txt \
                     -g ${HOME}/data/gene_coding.txt \
                     -s ${HOME}/TMPDATA/test_strains.txt \
                     -o ${HOME}/TMPDATA/test.SNPs.hb.txt

# statistical testing with trait data
bin/haplomap ghmap -p ${HOME}/data/test_traits.txt \
                   -b ${HOME}/TMPDATA/test.SNPs.hb.txt \
                   -o ${HOME}/TMPDATA/test.final.output.txt
```


## Input
1. eblocks:
    - Strain file (-s): column1 -> abbrev, column2 -> fullname
    - Allele file (-a): NIEHS compact format
    - Gene Annotation (-g): 
       - format: <SNP_{chr}_{postion}>  <gene_name>  < SNP_cateogry> 

2. ghmap:
    - Trait file (-p):  column1 -> abbrev, column2 -> value
    - eblock output (-o): haplotype blocks

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
|7 |popFDR       |   Pillai’s Trace FDR |
|8 |popYes       | Pillai’s Trace FDR Rejection | 
|9 |Chrom        | Chromosome idx      |
|10 |ChrBeg      | HaploBlock begin idx|
|11 |ChrEnd      | HaploBlock end idx  |
|12 |GeneExprMap | Gene expression Map |

  * block-oriented result file

BlockID | (IGNORED) | blockStart | blockSize | ChrIdx | ChrStart | ChrEnd | Pattern | Fstat/Pval | CondingMap ...

## TODO
[X] MPD trait data curation  
[-] Unified commandline interface 

## Changelog
v0.1
* ghmap: population structure testing (Pillai’s Trace)
* ghmap: multiple hypothesis testing correction. Benjamini Hochberg procedure
* ghmap: support raw animaldata (individual data) input, which will give more statistical power..
* Refactor source code, Cross-platform support (Linux, MacOS)
* Refactor to modern C++11

v0.0
* the first version of haplomap