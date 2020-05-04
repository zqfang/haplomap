import os
from snakemake.shell import shell
############################# Required ###################################
# set output directory 
WORKSPACE = "/data/bases/shared/haplomap/results_mpd20200422"
workdir: WORKSPACE

# MPD trait ids 
#TRAIT_IDS = os.path.join(WORKSPACE, "test_ids.txt")
TRAIT_IDS = "/data/bases/shared/haplomap/new_test_ids.txt"
# ghmap input
TRAIT_DATA =  "/data/bases/shared/haplomap/strainmeans_old_byGender.csv"

# eblock input
STRAIN_ANNO = "/data/bases/shared/haplomap/PELTZ_20180101/Strains_20180101.csv"
SNPS_DIR = "/data/bases/shared/haplomap/PELTZ_20180101/SNPS"
GENE_ANNO = "/data/bases/shared/haplomap/PELTZ_20180101/gene_coding.txt"

# VCF input pattern (from samtools or GATK)
VCFs = "/data/bases/fangzq/VCFs/combined.chr{i}.snp.filter.vcf"

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
    input: HBCGM,#SNPDB

# rule vcf2strains:
#     input:  
#         # vcf = "VCFs/combined.{chrom}.raw.vcf", 
#         vcf = VCFs
#     output: 
#         temp("SNPs/chr{i}.strains.temp")
#     shell:
#         # NOTE: '\t' is default delim for cut
#         "head -n 1000 {input.vcf} | grep '^#CHROM' | "
#         "cut -f10-  > {output}"  
    
# rule vcf2niehs:
#     input:  
#         # vcf = "VCFs/combined.chr{i}.raw.vcf", 
#         vcf = VCFs,
#         strains = "SNPs/chr{i}.strains.temp"
#     output: 
#         protected("SNPs/chr{i}.txt")
#     params:
#         outdir= "SNPs",
#         chrom="{i}",
#         qual_samtools=50, 
#         heterzygote_cutoff = 20
#     script:
#         "scripts/vcf2NIEHS.py"

# rule pheno:
#     input: TRAIT_DATA
#     ouput: os.path.join(OUTPUT_DIR, "mpd.ids.txt")
#     shell: 
#         "cut -d, -f1 {input} | uniq | sed '1d' > {output.txt}"
rule traits: 
    output: temp(expand("MPD_{ids}/strain.{ids}.temp", ids=IDS))
    run:
        for out in output:
            shell("touch %s"%out)

rule strain2trait:
    input: 
        trait = TRAIT_DATA,
        strain = STRAIN_ANNO,
        ids = "MPD_{ids}/strain.{ids}.temp"
    output: 
        "MPD_{ids}/strain.{ids}.txt",
        "MPD_{ids}/trait.{ids}.txt",
    params:
        outdir = WORKSPACE,
        traitid = "{ids}"
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
        "{params.smkdir}/haplomap/build/bin/eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} -o {output.hb} "
        "-v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "MPD_{ids}/chr{i}.hblocks.txt",
        trait = "MPD_{ids}/trait.{ids}.txt"
    output: "MPD_{ids}/chr{i}.results.txt"
    params:
        smkdir = SNAKEMAKE_DIR,
        cat = "MPD_{ids}/trait.{ids}.categorical"
    log: "logs/MPD_{ids}.chr{i}.ghmap.log"
    run:
        categorical = "-c" if os.path.exists(params.cat) else ''
        cmd = "{params.smkdir}/haplomap/build/bin/ghmap %s "%categorical +\
              "-p {input.trait} -b {input.hb} -o {output} " +\
              "-n MPD_{wildcards.ids} -v > {log}"
        shell(cmd)
