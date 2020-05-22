
# HBCGM
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)


## Dependency 

Ubuntu 18.04.2 LTS
* GSL
* GCC >= 4.8

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


## Usage  
See more detail in ``haplomap`` subfolder: [Install](haplomap/README.md)

### 1. Prepare MPD trait id file. Each id per row. e.g.
```
1501-f
26720-m
26720-f
...
```

### 2. Edit the required file path in `haplomap.smk`, including

```python
# working directory
WORKSPACE = "/data/bases/fangzq/20200429"

# MPD trait id file
TRAIT_IDS = "/data/bases/shared/haplomap/new_test_ids.txt"

# MPD trait database 
TRAIT_DATA =  "/data/bases/shared/haplomap/strainmeans_old_byGender.csv"

# strain metadata for SNP database
STRAIN_ANNO = "/data/bases/shared/haplomap/PELTZ_20180101/Strains_20180101.csv"

# SNP database
SNPS_DIR = "/data/bases/shared/haplomap/PELTZ_20180101/SNPS"
# genen annotation input 
GENE_ANNO = "/data/bases/shared/haplomap/PELTZ_20180101/gene_coding.txt"
```

### 3. run
```shell
# modify the file path in haplomap and run with 24 cores
snakemake -s haplomap.smk -k -p -j 24   
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

3. Get SNP database
  - select snakemake pipelines in workflow folder.
  - edit the input and output files, then run  
    e.g.
    ```shell
    # modify the file path in haplomap and run with 24 cores
    snakemake -s workflows/bcftools.call.smk -k -p -j 24   
    ```

## Output

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