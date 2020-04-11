import os
############################# Required ###################################
# set output directory 
WORKSPACE = "/data/bases/fangzq"
workdir: WORKSPACE

# VCF input pattern (from samtools or GATK)
VCFs = "/data/bases/fangzq/VCFs/combined.chr{i}.snp.filter.vcf"

# MPD trait ids 
TRAIT_IDS = "/data/bases/shared/haplomap/test_ids.txt"

# eblock input
STRAIN_ANNO = "/data/bases/shared/haplomap/PELTZ_20190301/Strains_20190301.csv"
SNPS_DIR = "/data/bases/shared/haplomap/PELTZ_20190301/SNPS"
GENE_ANNO = "/data/bases/shared/haplomap/PELTZ_20190301/gene_coding.txt"

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
SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOMES)
HBCGM =  expand("MPD_{ids}/chr{i}.results.txt", ids=IDS, i=CHROMOSOMES)
HBLOCKS = expand("MPD_{ids}/chr{i}.hblocks.txt", ids=IDS, i=CHROMOSOMES)

rule target:
    input: SNPDB, #HBCGM

rule vcf2strains:
    input:  
        # vcf = "VCFs/combined.{chrom}.raw.vcf", 
        vcf = VCFs
    output: 
        temp("SNPs/chr{i}.strains.txt")
    shell:
        # NOTE: '\t' is default delim for cut
        "head -n 1000 {input.vcf} | grep '^#CHROM' | "
        "cut -f10-  > {output}"  
    
rule vcf2niehs:
    input:  
        # vcf = "VCFs/combined.chr{i}.raw.vcf", 
        vcf = VCFs,
        strains = "SNPs/chr{i}.strains.txt"
    output: 
        protected("SNPs/chr{i}.txt")
    params:
        outdir= "SNPs",
        chrom="{i}",
        qual_samtools=50, 
        heterzygote_cutoff = 20
    script:
        "scripts/vcf2NIEHS.py"

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
        snps = os.path.join(SNPS_DIR, "chr{i}.txt"),
        gene_anno = GENE_ANNO,
        strains = "MPD_{ids}/strain.{ids}.txt",
    output: 
        hb = protected("MPD_{ids}/chr{i}.hblocks.txt"),
        snphb = temp("MPD_{ids}/chr{i}.snp.hblocks.txt")
    params:
        smkdir = SNAKEMAKE_DIR
    log: "logs/MPD_{ids}.chr{i}.eblocks.log"
    shell:
        "{params.smkdir}/build/bin/eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} -o {output.hb} "
        "-v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "MPD_{ids}/chr{i}.hblocks.txt",
        trait = "MPD_{ids}/trait.{ids}.txt"
    output: "MPD_{ids}/chr{i}.results.txt"
    params:
        smkdir = SNAKEMAKE_DIR
    log: "logs/MPD_{ids}.chr{i}.ghmap.log"
    shell:
        "{params.smkdir}/build/bin/ghmap -p {input.trait} -b {input.hb} "
        "-o {output} -v > {log}"
