
# Gene annotation file generation

Gene annotation file includes which genes a haplotype block is associated with and their semantic effects. This needs to be generated dynamically in the block formation step since the semantic effect varies depending on the strains selected. 

## Prerequisites 

1. Download ANNOVAR main package: https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/ (registration required)
2. Download appropriate databases from UCSC genome browser (as of 08/04/2020 ANNOVAR did not have a downloadble package for mm10 GENCODE gene set)
2.a. From https://genome.ucsc.edu/ click "Table Browser" under "Tools" tab
2.b. Select mm10 as assembly and "Genes and Gene Predictions" in group track "GENCODE VM23"
2.c. Download 2 tables: "knownGene" save as "mm10\_knownGene.txt" and "kgXref" save as "mm10\_kgXref.txt" (output format all fields from selected table)
2.d. Download mRNA sequence file by selecting "knownGene" table then selecting sequence as output format. "get output" will prompt the page to the next selection. Select "mRNA" and "submit" then "get sequence" on the next page with default parameters. Save as "mm10\_knownGeneMrna.fa"  
2.e. store all of these in 1 directory. This directory will be the first argument in the next command.
3. run "conda install tqdm" to install tqdm python package

## ANNOVAR and merging the output to Amino Acid changes by strains file

To run gene annotation by ANNOVAR use the following command. VCFs must be stored in 1 directory with filename extension "vcf" and vcf must start with the chromosome name (or identifying name). e.g. chr1.vcf.
After annotated files are generated, these files will be merged into one gene annotation file, then split amino acid changes by strains. [WARNING] this step will take a long time. 

```bash
sh annovar_pipe_hbcgm.sh /path/to/files/from/prerequisites/step/ /path/to/vcf /path/to/ANNOVAR /path/to/SNPs/file/for/HBCGM
```

Final output will be saved in "AA\_by\_strains\_\*.pkl" split by chromosome. This file along with the mm10\_kgXref.txt and mm10\_knownGene.txt will be used as inputs for HBCGM (update the path in config.yaml file)



#example


