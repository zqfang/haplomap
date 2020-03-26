import os
############################# Required ###################################
# set output directory 
WORKSPACE = "/data/bases/shared/haplomap/test_results_20200324"
workdir: WORKSPACE

# MPD trait ids 
TRAIT_IDS = "/data/bases/shared/haplomap/test_ids2.txt"

# eblock input
STRAIN_ANNO = "/data/bases/shared/haplomap/Strains_20180101.csv"
SNPS_DIR = "/data/bases/shared/haplomap/SNPs_20180101"
GENE_ANNO = "/data/bases/shared/haplomap/gene_coding.txt"

# ghmap input
TRAIT_DATA =  "/data/bases/shared/haplomap/test_strainmeans_20200324.csv"

############################################################################
# snakefile dir
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

## trait ids
with open(TRAIT_IDS, 'r') as t:
    IDS = t.read().strip().split()

CHROMOSOMES = [str(i) for i in range (1, 20)] + ['X'] # NO 'Y'
# output files
HBCGM =  expand("MPD_{ids}/chr{chrm}.HBCGM.txt", ids=IDS, chrm=CHROMOSOMES)
HBLOCKS = expand("MPD_{ids}/{chrm}.hblocks.txt", ids=IDS, chrm=CHROMOSOMES)

rule target:
    input: HBCGM

# rule pheno:
#     input: TRAIT_DATA
#     ouput: os.path.join(OUTPUT_DIR, "mpd.ids.txt")
#     shell: 
#         "cut -d, -f1 {input} | uniq | sed '1d' > {output.txt}"

rule strain2trait:
    input: 
        trait = TRAIT_DATA,
        strain = STRAIN_ANNO,
        ids = TRAIT_IDS
    output: 
        #trait = directory(OUTPUT_DIR)
        expand("MPD_{ids}/strain.{ids}.txt", ids=IDS),
        expand("MPD_{ids}/trait.{ids}.txt", ids=IDS)
    params:
        outdir = WORKSPACE
    script:
        "scripts/strain2traits.py"

# find haplotypes
rule eblocks:
    input: 
        snps = os.path.join(SNPS_DIR, "{chrm}.txt"),
        gene_anno = GENE_ANNO,
        strains = "MPD_{ids}/strain.{ids}.txt",
    output: 
        hb = protected("MPD_{ids}/{chrm}.hblocks.txt"),
        snphb = temp("MPD_{ids}/{chrm}.snp.hblocks.txt")
    params:
        smkdir = SNAKEMAKE_DIR
    log: "logs/MPD_{ids}.{chrm}.eblocks.log"
    shell:
        "{params.smkdir}/build/bin/eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} -o {output.hb} "
        "-v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "MPD_{ids}/{chrm}.hblocks.txt",
        trait = "MPD_{ids}/trait.{ids}.txt"
    output: 
         "MPD_{ids}/{chrm}.results.txt"
    params:
        smkdir = SNAKEMAKE_DIR
    log: "logs/MPD_{ids}.{chrm}.ghmap.log"
    shell:
        "{params.smkdir}/build/bin/ghmap -p {input.trait} -b {input.hb} "
        "-o {output} -v > {log}"
