
# HBCGM
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)

## Installation

1. Install GSL and export the lib path  
e.g.
```bash
./configure --prefix=${HOME}/program/gsl
make && make install
```

2. edit `CMakeLists.txt`, set GSL header and lib path, 

```cmake
set(GSL_INCLUDE /path/to/gsl/include)
set(GSL_LIBS /path/to/gsl/lib)
```


3. build
```bash
mkdir build && cd build
cmake ..
make
```

## Dependency 

Ubuntu 18.04.2 LTS
* GSL 2.6
* GCC 7.4.0

## Usage
e.g.
```bash
# find haplotypes
eblocks -a ${HOME}/data/SNPS/chr18.txt \
        -g ${HOME}/data/gene_coding.txt \
        -s ${HOME}/TMPDATA/test_strains.txt \
        -p ${HOME}/TMPDATA/test.haploblocks.txt \
        -o ${HOME}/TMPDATA/test.SNPs.hb.txt

# statistical testing with trait data
ghmap -p ${HOME}/data/test_traits.txt \
      -b ${HOME}/TMPDATA/test.SNPs.hb.txt \
      -o ${HOME}/TMPDATA/test.final.output.txt
```

### Snakemake pipeline  
1. Prepare MPD trait id file. Each id per row
2. Edit the required file path in `haplomap.smk`
3. run
```shell
# modify the file path in haplomap and run with 24 cores
snakemake -s haplomap.smk -p -j 24  
```


### Input
1. eblocks:
    - Strain file (-s): column1 -> abbrev, column2 -> fullname
    - Allele file (-a): NIEHS compact format
    - Gene Annotation (-g): 
       - format: <SNP_{chr}_{postion}>  <gene_name>  < SNP_cateogry> 

2. ghmap:
    - Trait file (-p):  column1 -> abbrev, column2 -> value
    - eblock output (-o): haplotype blocks

### Output

1. ebloks:

|NO|Field | Explaination |
|--| ---- | ------------ |
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

|NO|Field | Explaination |
|--| ---- | ------------ |
|0 |GeneName   | Associated Gene     |
|1 |CodonFlag  | condonChange ? 1:0  |
|2 |Pattern    | Haplotype pattern   |
|3 |Score      | isCategorical ? Fstat : Pvalue |
|4 |Effect     | Genetic Effect (Effect Size)   |
|5 |Chrom      | Chromosome idx      |
|6 |ChrBeg     | HaploBlock begin idx|
|7 |ChrEnd     | HaploBlock end idx  |
|8 |GeneExprMap| Gene expression Map |

  * block-oriented result file

BlockID | (IGNORED) | blockStart | blockSize | ChrIdx | ChrStart | ChrEnd | Pattern | Fstat/Pval | CondingMap ...