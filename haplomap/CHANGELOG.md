# haplomap
Haplotype-based computational genetic mapping  

![Haplomap](https://github.com/zqfang/haplomap/workflows/Haplomap/badge.svg)



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
  - minor patches for structural variant input.
* pca: 
  - add a new sub-command 
  - could be used for getting genetic relationship 
* covnert:
  - For historical reasons, eblocks use NIEHS compact format as input. 
  - To use haplomap more friendly, we now do this for you. 
  - support structral variant input. 
* annotate:
  - convert ensemble-vep result for variant annotation input for eblocks.
 

v0.0
* The original version of haplomap