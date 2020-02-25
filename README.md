
# HBCGM
haplotype-based computational genetic mapping

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
ghmap -p ${HOME}/data/compact_gene_expr.txt \
      -b ${HOME}/TMPDATA/test.SNPs.hb.txt \
      -o ${HOME}/TMPDATA/test.final.ouput.txt
```

### Some Questions
1. If SNPs pattern has 'D' (deletion) in all input strains, drop this SNP 
2. Strain file: column1 -> name, column2 -> abbrev

