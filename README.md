
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
# you may need to add this line to your .bashrc 
export LD_LIBRARY_PATH="${HOME}/program/gsl/lib:$LD_LIBRARY_PATH"
```

2. edit `CMakeLists.txt`, set GSL header and lib path, 

```cmake
set(GSL_INCLUDE /path/to/gsl/include)
set(GSL_LIBS /path/to/gsl/lib)
```

3. build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```


## Usage  

### 1. Run haplomap standalone
See more detail in ``haplomap`` subfolder: [Run haplomap standalone](haplomap/README.md)

### 2. Use `snakemake` workflow to run
#### (1). Prepare MPD `measnum` id file. One id per row, suffixed with "-m" or "-f"(f: female, m: male)
```
26720-m
26720-f
9940-f
...


```

#### (2). Edit the `config.yaml` file path in `workflows` folder:

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

#### (3). run haplomap pipeline

Install `snakemake` first. You need `Miniconda` if conda is not installed

#### (3.1) create conda envs
```shell
conda create -n hbcgm -f environment.yaml
```

#### (3.2) run on a local computing node.

```shell
source activate hbcgm
# modify the file path in haplomap and run with 24 cores
snakemake -s workflows/haplomap.smk \
          --configfile workflows/config.yaml 
          -k -p -j 24   
```
### Run on the HPC, e.g. Stanford Sherlock 

e.g. Sherlock slurm
1. edit `slurm.submit.sh`, change file path to `HBCGM/workflows`
2. edit `workflows/slurm_config.yaml`, specify the resource you need.
3. submit
```
sbatch slurm.submit.sh
```


## Output
output explanation, see here: [Run haplomap standalone](haplomap/README.md)



## Contact

Email: 
- Zhuoqing Fang: fangzq@stanford.edu
- Gary Peltz: gpeltz@stanford.edu

