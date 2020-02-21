# HBCGM
haplotype-based computational genetic mapping

This is the **original version** of HBCGM

## Installation

1. Install GSL and export the lib path  
e.g.
```bash
./configure --prefix=${HOME}/program/gsl
make && make install
```
2. edit the Makefile -L /path/to/gsl/lib -I /path/to/gsl/include.  
e.g.
```bash
g++ $(CFLAGS) -c \
    -L${HOME}/program/gsl/program/gsl/lib -l gsl -l gslcblas \
    -I${HOME}/program/gsl/program/gsl/include -I./include \
    quantTraitMap.cpp
```
3. run the Makefile
```bash
cd haplomap
make eblocks
make ghmap
```

