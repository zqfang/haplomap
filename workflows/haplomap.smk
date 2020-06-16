import os, glob
############################# Required ###################################
# set output directory 
#configfile: "config.yaml"
workdir: config['HBCGM']['WORKSPACE']
HBCGM_BIN = config['HBCGM']['BIN']
# MPD trait ids 
TRAIT_IDS = config['HBCGM']['TRAIT_IDS']
# ghmap input
TRAIT_DATA =  config['HBCGM']['TRAIT_DATA']
# eblock input
STRAIN_ANNO = config['HBCGM']['STRAIN_ANNO']
SNPDB = config['HBCGM']['SNPS_DIR']
GENE_ANNO = config['HBCGM']['GENE_ANNO']
GENE_EXPRS = config['HBCGM']['GENE_EXPRS']
# open chromatin regions input
ATAC_PEAKS = glob.glob(config['HBCGM']['ATAC_PEAKS'])

############################################################################

## trait ids
with open(TRAIT_IDS, 'r') as t:
    IDS = t.read().strip().split()
    assert len(IDS) >= 1

CHROMOSOMES = [str(i) for i in range (1, 20)] + ['X'] # NO 'Y'
# output files
# SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOMES)
HBCGM =  expand("MPD_{ids}/chr{i}.results.txt", ids=IDS, i=CHROMOSOMES)
HBLOCKS = expand("MPD_{ids}/chr{i}.hblocks.txt", ids=IDS, i=CHROMOSOMES)
HBCGM_NONCODING = expand("MPD_{ids}/chr{i}.open_region.bed", ids=IDS, i=CHROMOSOMES)
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
        outdir = config['HBCGM']['WORKSPACE'],
        traitid = "{ids}",
        has_categorical = not config['HBCGM']['NUMERIC']
    script:
        "../scripts/strain2traits.py"

# find haplotypes
rule eblocks:
    input: 
        snps = os.path.join(SNPDB, "chr{i}.txt"),
        gene_anno = GENE_ANNO,
        strains = "MPD_{ids}/strain.{ids}.txt",
    output: 
        hb = protected("MPD_{ids}/chr{i}.hblocks.txt"),
        snphb = temp("MPD_{ids}/chr{i}.snp.hblocks.txt")
    params:
        bin = HBCGM_BIN,
    log: "logs/MPD_{ids}.chr{i}.eblocks.log"
    shell:
        "{params.bin}/eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} "
        "-o {output.hb} -v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "MPD_{ids}/chr{i}.hblocks.txt",
        trait = "MPD_{ids}/trait.{ids}.txt",
        gene_exprs = GENE_EXPRS,
    output: "MPD_{ids}/chr{i}.results.txt"
    params:
        bin = HBCGM_BIN,
        cat = "MPD_{ids}/trait.{ids}.categorical"
    log: "logs/MPD_{ids}.chr{i}.ghmap.log"
    run:
        categorical = "-c" if os.path.exists(params.cat) else ''
        cats = "catogorical" if os.path.exists(params.cat) else ''
        cmd = "{params.bin}/ghmap %s "%categorical +\
              "-e {input.gene_exprs} " +\
              "-p {input.trait} -b {input.hb} -o {output} " +\
              "-n MPD_{wildcards.ids}_%s -v > {log}"%cats
        shell(cmd)

rule open_chrom:
    input: 
        ghmap="MPD_{ids}/chr{i}.results.txt",
        peaks=ATAC_PEAKS,
    output:
        "MPD_{ids}/chr{i}.open_region.bed",
    params:
        peaks=" ".join(ATAC_PEAKS)
    shell:
        "sed '1,3d' {input.ghmap} | cut -f6-8 | "
        # add 'chr' and convert hblocks to 0-based coordinate #BEGIN {{FS = \"\\t\";OFS = \"\\t\" }};
        "awk -F'\\t' -v OFS='\\t' -v s=1 '{{print \"chr\"$1, $2-s, $3}}' | "
        "sort -k1,1 -k2,2n | "
        "bedtools intersect -sorted -a stdin -b {params.peaks} -wa -wb > {output}"

# about sort -k field1[,field2]
# To sort on the first field and then on the second: sort -k1,1 -k2,2
# The arguments field1 and field2 have the form m.n (m,n > 0) 
# and can be followed by one or more of the modifiers b, d, f, i, n, g, M and r,
# which correspond to the sort options. -n: numberic sort


