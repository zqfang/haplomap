import os, sys, json


WKDIR = "/home/fangzq/github/HBCGM/example/PeltzData"

workdir: WKDIR
# INPUT_FILES
MESH_DICT = {'10806':'D003337',
 '10807':'D020712',
 '10813':'D020712',
 '26721': 'D017948',
 '50243':'D002331',
 '9904': 'D008076',
 'Haloperidol': "D064420",
 'irinotecan': 'D000077146',
 'Cocaine':'D007858',
}

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
MPD2MeSH = "/data/bases/shared/haplomap/PELTZ_20210609/mpd2mesh.json"
# with open(MPD2MeSH, 'r') j:
#     mpd2mesh = json.load(j)

### OUTPUT

CHROMOSOMES = [str(i) for i in range (1, 20)] + ['X'] # NO 'Y'
# output files
# SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOMES)
HBCGM =  expand("MPD_{ids}/chr{i}.results.txt", ids = list(MESH_DICT.keys()), i=CHROMOSOMES)
MESH = expand("MPD_{ids}.results.mesh.txt", ids=list(MESH_DICT.keys()))

rule target:
    input: HBCGM, MESH


rule annotateSNPs:
    input:
        strains = "Peltz_{ids}.trait.txt",
        snps = os.path.join(SNPDB_DIR, "chr{i}.txt"), 
        annodb = os.path.join(ANNOVAR, "chr{i}.AA_by_strains.pkl"),
        kgxref = KNOWNGENE_META, 
        knowngene= KNOWNGENE,
    output:
        hgnc = temp("MPD_{ids}/hblocks/chr{i}.genename.txt"),
        ensemble = temp("MPD_{ids}/hblocks/chr{i}.geneid.txt"),
    script:
        "../scripts/annotateSNPs.py"



# find haplotypes
rule eblocks:
    input: 
        snps = os.path.join(SNPDB_DIR, "chr{i}.txt"),
        gene_anno = "MPD_{ids}/hblocks/chr{i}.genename.txt",# "MPD_{ids}/gene_coding.txt",#
        strains = "Peltz_{ids}.trait.txt",
    output: 
        hb = "MPD_{ids}/hblocks/chr{i}.hblocks.txt",
        snphb ="MPD_{ids}/hblocks/chr{i}.snp.hblocks.txt"
    params:
        bin = HBCGM_BIN,
    log: "logs/MPD_{ids}.chr{i}.eblocks.log"
    shell:
        "{params.bin}/haplomap eblocks -a {input.snps} -g {input.gene_anno} "
        "-s {input.strains} -p {output.snphb} "
        "-o {output.hb} -v > {log}"

# statistical testing with trait data       
rule ghmap:
    input: 
        hb = "MPD_{ids}/hblocks/chr{i}.hblocks.txt",
        trait = "Peltz_{ids}.trait.txt",
        gene_exprs = GENE_EXPRS,
        rel = GENETIC_REL,
    output: "MPD_{ids}/chr{i}.results.txt"
    params:
        bin = HBCGM_BIN,
        cat = "MPD_{ids}/trait.{ids}.categorical"
    log: "logs/MPD_{ids}.chr{i}.ghmap.log"
    run:
        categorical = "-c" if os.path.exists(params.cat) else ''
        cats = "_catogorical" if os.path.exists(params.cat) else ''
        cmd = "{params.bin}/haplomap ghmap %s "%categorical +\
              "-e {input.gene_exprs} -r {input.rel} " +\
              "-p {input.trait} -b {input.hb} -o {output} " +\
              "-n MPD_{wildcards.ids}%s -a -v > {log}"%cats
        shell(cmd)

rule mesh:
    input: ["MPD_{ids}/chr%s.results.txt"%c for c in CHROMOSOMES]
    output: "MPD_{ids}.results.mesh.txt"
    params: 
        mesh = lambda wildcards: MESH_DICT[wildcards.ids]
    threads: 1
    shell:
        "/home/fangzq/miniconda/envs/fastai/bin/python "
        "/home/fangzq/github/InpherGNN/GNNphrank/predict.py "
        "--bundle /data/bases/fangzq/Pubmed/bundle " 
        "--hbcgm_result_dir MPD_{wildcards.ids} "
        "--mesh_terms {params.mesh} --num_cpus {threads} "