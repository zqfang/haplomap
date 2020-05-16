import os, glob
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

# open chromatin regions input
ATAC_PEAKS = glob.glob("/data/bases/fangzq/MouseEpigenomeAtlas/beds/*.blacklist_removed.broadPeak")

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
HBCGM_NONCODING = expand("MPD_{ids}/chr{i}.regulatory_region.results.txt", ids=IDS, i=CHROMOSOMES)
rule target:
    input: HBCGM, HBCGM_NONCODING


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

rule open_chrom:
    input: 
        ghmap="MPD_{ids}/chr{i}.results.txt",
        peaks=ATAC_PEAKS,
    output:
        "MPD_{ids}/chr{i}.regulatory_region.results.txt",
    params:
        peaks=" ".join(ATAC_PEAKS)
    shell:
        "sed '1,3d' {input.ghmap} | cut -f6-8 | "
        "sed -e 's/^/chr/' | sort -k1,1 -k2,2n | "
        "bedtools intersect -sorted -a stdin -b {params.peaks} -wa -wb > {output}"

# about sort -k field1[,field2]
# To sort on the first field and then on the second: sort -k1,1 -k2,2
# The arguments field1 and field2 have the form m.n (m,n > 0) 
# and can be followed by one or more of the modifiers b, d, f, i, n, g, M and r,
# which correspond to the sort options. -n: numberic sort


