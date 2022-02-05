
# Haplomap (a.k.a HBCGM )
Haplotype-based computational genetic mapping  
![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)

This is a brand-new version of HBCGM, named Haplomap. 

## Dependency 
Works both on `Linux`and `MacOS`

Haplomap:
* CMake
* GCC >= 4.8
* clang >= 11.0.3 (only tested with 11.x version)
* C++11
* GSL

For Variant Calling, you need:
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

### Run haplomap standalone
See more detail in ``haplomap`` subfolder: [Run haplomap standalone](haplomap/README.md)

### Use `snakemake` workflow to run Mouse Phenome Database (MPD) datasets

[Mouse Phenome Database](https://phenome.jax.org/) have > 10K datasets. Try to configure the files below to run
#### 1. Prepare MPD `measnum` id file. One id per row, suffixed with "-m" or "-f"(f: female, m: male)
```
26720-m
26720-f
9940-f
...
```

#### 2. Edit the `config.yaml` file path in `workflows` folder:

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
    # set to true select individual animal data. Default: use strain means.   
    USE_RAWDATA: false 
    # strains metadata: map strain abbrev to full name, jax ids, ect. 
    # see docs folder to view example
    STRAIN_ANNO: "/data/bases/shared/haplomap/PELTZ_20210609/strains.metadata.csv"

    # eblock input
    # `haplomap convert` output after variant calling step 
    SNPS_DIR: "/data/bases/shared/haplomap/PELTZ_20210609/SNPs"
    # Ensembl-vep output after variant calling step
    VEP_DIR: "/data/bases/shared/haplomap/PELTZ_20210609/VEP"

    ## Optional files
    # genetic relation file from PLink output
    GENETIC_REL: "/data/bases/shared/haplomap/PELTZ_20210609/mouse54_grm.rel"
    # gene expression file 
    GENE_EXPRS: "/data/bases/shared/haplomap/PELTZ_20210609/mus.compact.exprs.txt"
```

#### 3. run haplomap pipeline

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
#### 3.3 Run on the HPC, e.g. Stanford Sherlock 

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

## Copyright and License Information
Copyright (C) 2019-2022 Stanford University, Zhuoqing Fang and Gary Peltz.

Authors: Zhuoqing Fang and Gary Peltz

This program is licensed with commercial restriction use license. Please see the attached LICENSE file for details.