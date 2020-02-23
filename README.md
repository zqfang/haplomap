
# HBCGM
haplotype-based computational genetic mapping

## Installation

1. Install GSL and export the lib path  
e.g.
```bash
./configure --prefix=${HOME}/program/gsl
make && make install
```

2. edit `CMakeLists.txt`, set GSL header and lib path, then
```bash
mkdir build && cd build
cmake ..
make
```

3. Alternatively, edit the Makefile, set -L /path/to/gsl/lib -I /path/to/gsl/include.  
e.g.
```bash
g++ $(CFLAGS) -c \
    -L${HOME}/program/gsl/program/gsl/lib -l gsl -l gslcblas \
    -I${HOME}/program/gsl/program/gsl/include -I./include \
    quantTraitMap.cpp
```
the run the Makefile
```bash
cd haplomap
make eblocks
make ghmap
```

## Dependency 

Ubuntu 18.04.2 LTS
* GSL 2.6
* GCC 7.4.0


### Some tips
1. If SNPs pattern has 'D' (deletion) in all input strains, drop this SNP 

