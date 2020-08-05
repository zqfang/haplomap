
# HBCGM
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)


## Dependency 
Works both on `Linux`and `MacOS`

Haplomap:
* CMake
* GCC >= 4.8
* clang >= 11.0.3 (only tested with 11.x version)
* C++11
* GSL

Germline Variant Calling
* GATK 4.x
* SAMtools
* BCFtools
* BEDtools
* BWA

Running pipeline
* Snakemake


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

### 1. Prepare MPD `measnum` id file. Each id per row. Suffix for sex(f,m) permitted. e.g.
```
1501
03504
26720-m
26720-f
...
```

### 2. Edit the `config.yaml` file path in `workflows` folder:

only edit `HBCGM` section.
```yaml
HBCGM:
    # working directory
    WORKSPACE: "/data/bases/fangzq/MPD/results_drug_diet"
    # path to haplomap
    BIN: "/home/fangzq/github/HBCGM/build/bin"
    
    # ghmap input
    # MPD trait ids 
    TRAIT_IDS: "/data/bases/fangzq/MPD/drug-diet.ids.txt"
    # set to true if input individual animal data. Default: use strain means.   
    USE_RAWDATA: false 
    # given file path, use input data instead of using MPD API to get data.
    TRAIT_DATA:  "" #"/data/bases/shared/haplomap/AHresponse_strainmeans2.txt"
    # genetic relation file from PLink output
    GENETIC_REL: "/data/bases/shared/haplomap/PELTZ_20200429/mouse54_grm.rel"
    GENETIC_REL_ID: "/data/bases/shared/haplomap/PELTZ_20200429/mouse54_grm.rel.id"
    # gene expression file
    GENE_EXPRS: "/data/bases/shared/haplomap/PELTZ_20200429/mus.compact.exprs.txt"
    # open chromatin regions to annotate haploblocks.
    ATAC_PEAKS: "/data/bases/fangzq/MouseEpigenomeAtlas/beds/*.blacklist_removed.broadPeak"

    # eblock input
    # strains metadata. 
    STRAIN_ANNO: "/data/bases/shared/haplomap/PELTZ_20200429/strains.metadata.csv"
    # path to SNP database
    SNPS_DIR: "/data/bases/shared/haplomap/PELTZ_20200429/SNPs"
    # SNP annotations for all genes
    ANNOVAR: "/data/bases/shared/haplomap/PELTZ_20200429/AA_by_strains.pkl" 
    KNOWNGENE_META: "/data/bases/shared/haplomap/PELTZ_20200429/mm10_kgXref.txt" 
    KNOWNGENE: "/data/bases/shared/haplomap/PELTZ_20200429/mm10_knownGene.txt" 
```

### 3. run haplomap pipeline

Install `snakemake` first. You need `Miniconda` if conda is not installed

#### 3.1 create conda envs
```shell
conda create -n hbcgm -f environment.yaml
```

#### 3.2 run on a local computing node.

```shell
source activate hbcgm
# modify the file path in haplomap and run with 24 cores
snakemake -s workflows/haplomap.smk \
          --configfile workflows/config.yaml 
          -k -p -j 24   
```
#### 2.2 run on the HPC 

e.g. Sherlock slurm
1. edit `slurm.submit.sh`, change file path to `HBCGM/workflows`
2. edit `workflows/slurm_config.yaml`, specify the resource you need.
3. submit
```
sbatch slurm.submit.sh
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

| NO |Field | Explaination |
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