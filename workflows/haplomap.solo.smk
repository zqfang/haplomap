import os, sys


WKDIR = "/data/bases/fangzq/Pubmed/test_cases/TEST"

workdir: WKDIR
# INPUT_FILES

IDS = ['000']
STRAINS = "/data/bases/fangzq/Pubmed/test_cases/TEST/strain.txt"
STRAINS_VALUE = "/data/bases/fangzq/Pubmed/test_cases/TEST/strain.values.txt"
# genetic relation file from PLink output
GENETIC_REL =  "/data/bases/shared/haplomap/PELTZ_20210609/mouse54_grm.rel"
# gene expression file
GENE_EXPRS = "/data/bases/shared/haplomap/PELTZ_20210609/mus.compact.exprs.txt"
######### eblock input ######### 
HBCGM_BIN = "/home/fangzq/github/HBCGM/build/bin"
# strains metadata. 
STRAIN_ANNO = "/data/bases/shared/haplomap/PELTZ_20210609/strains.metadata.csv"
# path to SNP database
SNPDB_DIR =  "/data/bases/shared/haplomap/PELTZ_20210609/SNPs"
# PATH to SNP annotations for all genes
ANNOVAR = "/data/bases/shared/haplomap/PELTZ_20210609/SNP_Annotation" 
# snp, geneid,genename mapping
KNOWNGENE_META = "/data/bases/shared/haplomap/PELTZ_20210609/SNP_Annotation/mm10_kgXref.txt" 
KNOWNGENE = "/data/bases/shared/haplomap/PELTZ_20210609/SNP_Annotation/mm10_knownGene.txt" 

### OUTPUT

CHROMOSOMES = [str(i) for i in range (1, 20)] + ['X'] # NO 'Y'
# output files
# SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOMES)
HBCGM =  expand("RUN_{ids}/chr{i}.results.txt", ids = IDS, i=CHROMOSOMES)
HBLOCKS = expand("RUN_{ids}/chr{i}.hblocks.txt", ids = IDS, i=CHROMOSOMES)



rule target:
    input: HBCGM


rule annotateSNPs:
    input:
        strains = "RUN_{ids}/strain.{ids}.txt",
        snps = os.path.join(SNPDB_DIR, "chr{i}.txt"), 
        annodb = os.path.join(ANNOVAR, "chr{i}.AA_by_strains.pkl"),
        kgxref = KNOWNGENE_META, 
        knowngene= KNOWNGENE,
    output:
        hgnc = "RUN_{ids}/annot/chr{i}.genename.txt",
        ensemble = "RUN_{ids}/annot/chr{i}.geneid.txt",
    script:
        "../scripts/annotateSNPs.py"



# find haplotypes
rule eblocks:
    input: 
        snps = os.path.join(SNPDB_DIR, "chr{i}.txt"),
        gene_anno = "RUN_{ids}/annot/chr{i}.genename.txt",# "RUN_{ids}/gene_coding.txt",#
        strains = "RUN_{ids}/strain.{ids}.txt",
    output: 
        hb = protected("RUN_{ids}/haploblock/chr{i}.hblocks.txt"),
        snphb ="RUN_{ids}/haploblock/chr{i}.snp.hblocks.txt"
    params:
        bin = HBCGM_BIN,
    log: "logs/RUN_{ids}.chr{i}.eblocks.log"
    shell:
        "{params.bin}/haplomap eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} "
        "-o {output.hb} -v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "RUN_{ids}/haploblock/chr{i}.hblocks.txt",
        trait = "RUN_{ids}/trait.{ids}.txt",
        gene_exprs = GENE_EXPRS,
        rel = GENETIC_REL,
    output: "RUN_{ids}/chr{i}.results.txt"
    params:
        bin = HBCGM_BIN,
        cat = "RUN_{ids}/trait.{ids}.categorical"
    log: "logs/RUN_{ids}.chr{i}.ghmap.log"
    run:
        categorical = "-c" if os.path.exists(params.cat) else ''
        cats = "catogorical" if os.path.exists(params.cat) else ''
        cmd = "{params.bin}/haplomap ghmap %s "%categorical +\
              "-e {input.gene_exprs} -r {input.rel} " +\
              "-p {input.trait} -b {input.hb} -o {output} " +\
              "-n RUN_{wildcards.ids}_%s -v > {log}"%cats
        shell(cmd)