#!bin/sh

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