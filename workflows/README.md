# Snakemake pipeline for SNV calling


### Two pipeline developed:
1. bcftools call 
  - prefer pipeline for inbred mouse and HBCGM input
  - more accuary for inbred mouse


2. GATK best practise
  - designed for human genetics 
  - have to play with ``VQSR`` or ``hardfilering`` parameters if use non-human data


**Caution !**: Both pipeline takes a long time to run.


### Why not GATK

One of my colleague who studies mouse genetics, said, 

> I tried the haplotype caller from GATK. But it seems that the haplotype caller is designed for heterogeneous genome like human than for mice. Therefore, the result coming out of HC is worse than samtools, as I manually inspected a few regions that HC calls didn't make sense.

> In addition, in one of their mouse genomic paper that we reviewed, they even skipped the second recalibration step. We asked them why and they said it was because of the same reason: good for human but not that good for the homogeneous inbred mouse.