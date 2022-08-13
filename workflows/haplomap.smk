import os, glob, json
import pandas as pd
############################# Required ###################################
# set output directory 
#configfile: "config.yaml"
workdir: config['HBCGM']['WORKSPACE']
HBCGM_BIN = config['HBCGM']['BIN']
# MPD trait ids
if 'TRAIT_IDS' in config:
    TRAIT_IDS = config['TRAIT_IDS'] # input by --config TRAIT_IDS="ids.txt"
else:
    TRAIT_IDS = config['HBCGM']['TRAIT_IDS']

# ghmap input
TRAIT_DATA =  config['HBCGM']['TRAIT_DATA']
GENETIC_REL = config['HBCGM']['GENETIC_REL']
MPD2MeSH = config['HBCGM']['MPD2MeSH'] # json file
# eblock input
STRAIN_ANNO = config['HBCGM']['STRAIN_ANNO']
SNPDB = config['HBCGM']['VAR_DIR']
VCF_DIR = config['HBCGM']['VCF_DIR']
VEP_DIR = config['HBCGM']['VEP_DIR']
# ANNOVAR = config['HBCGM']['ANNOVAR'] 
# KNOWNGENE_META = config['HBCGM']['KNOWNGENE_META']
# KNOWNGENE = config['HBCGM']['KNOWNGENE']
GENE_EXPRS = config['HBCGM']['GENE_EXPRS']
# open chromatin regions input
ATAC_PEAKS = glob.glob(config['HBCGM']['ATAC_PEAKS'])

############################################################################

## trait ids
with open(TRAIT_IDS, 'r') as t:
    IDS_ = t.read().strip().split()
    assert len(IDS_) >= 1

# filter id in MeSH
with open(MPD2MeSH, 'r') as j:
    MESH_DICT = json.load(j)
IDS = [i for i in IDS_ if i.split("-")[0] in MESH_DICT ]
# filter id that already run
pat = config['HBCGM']['WORKSPACE'] + "MPD_{ids}_snp.results.txt"
IDS = [mnum for mnum in IDS if not os.path.exists(pat.format(ids = mnum))]

CHROMSOME = [str(i) for i in range (1, 20)] + ['X'] # NO 'Y'
# output files
# SNPDB = expand("SNPs/chr{i}.txt", i=CHROMOSOMES)
HBCGM =  expand("MPD_{ids}/chr{i}.snp.results.txt", ids=IDS, i=CHROMSOME)
HBCGM_NONCODING = expand("MPD_{ids}/chr{i}.open_region.bed", ids=IDS, i=CHROMSOME)
HBCGM_MESH = expand("MPD_{ids}_snp.results.mesh.txt", ids=IDS)
# rules that not work in a new node
#localrules: target, traits, strain2trait  


rule target:
    input: HBCGM_MESH, #HBCGM_NONCODING

########################### Prepare Phenotypic DATA (from MPD API) ############################
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
        strain = STRAIN_ANNO,
        ids = "MPD_{ids}/strain.{ids}.temp"
    output: 
        "MPD_{ids}/trait.{ids}.txt",
    params:
        trait = TRAIT_DATA,
        outdir = config['HBCGM']['WORKSPACE'],
        traitid = "{ids}",
        rawdata = config['HBCGM']['USE_RAWDATA']
    script:
        "../scripts/strain2traits.py"


############################ Convert VCF to niehs, tped, tfam  ########################################

rule snp2NIEHS:
    input:  
        vcf = os.path.join(VCF_DIR, "chr{i}.vcf.gz"),
    output: 
        niehs = protected(os.path.join(SNPDB, "chr{i}.snp.txt")),
        tped = temp(os.path.join(SNPDB, "chr{i}.tped")),
        tfam = temp(os.path.join(SNPDB, "chr{i}.tfam")),
    params:
        qual = config['BCFTOOLS']['qual'], 
        het = config['BCFTOOLS']['phred_likelihood_diff'],
        ad = config['BCFTOOLS']['allele_depth'],
        ratio = config['BCFTOOLS']['allele_mindepth_ratio'],
        mq = config['BCFTOOLS']['mapping_quality'],
        sb = config['BCFTOOLS']['strand_bias_pvalue'], 
        BIN = config['HBCGM']['BIN']# path to haplomap binary
    log: "logs/chr{i}.snp2niehs.log"
    shell:
        "bcftools view -f .,PASS -v snps {input.vcf} | "
        "{params.BIN}/haplomap convert -o {output.niehs} -a {params.ad} -r {params.ratio} "
        "-q {params.qual} -d {params.het} -m {params.mq} -b {params.sb} "
        "--plink -v > {log}"

################## Generate Relationship Matrix #################################################
rule toPLinkBed:
    input:
        os.path.join(SNPDB, "chr{i}.tped"),
        os.path.join(SNPDB, "chr{i}.tfam"),
    output:
        temp(os.path.join(SNPDB, "chr{i}.bed")),
        temp(os.path.join(SNPDB, "chr{i}.bim")),
        temp(os.path.join(SNPDB, "chr{i}.fam")),
    params:
        prefix = os.path.join(SNPDB, "chr{i}")
    shell:
        "plink --tfile {params.prefix} --make-bed --out {params.prefix}"
        

rule merge_list:
    output: temp(os.path.join(SNPDB, "mergelist.txt"))
    params:
        chrom = CHROMSOME,
        prefix = os.path.join(SNPDB, "chr")
    run:
        path = params.prefix
        with open(output[0], 'w') as out:
            for i in params.chrom:
                if str(i) == "1": continue 
                outline = f"{path}{i}.bed {path}{i}.bim {path}{i}.fam\n"
                out.write(outline)

rule merge_bed:
    input: 
        merge = os.path.join(SNPDB, "mergelist.txt"),
        beds = expand(os.path.join(SNPDB, "chr{i}.bed"), i=CHROMSOME),
        bims = expand(os.path.join(SNPDB, "chr{i}.bim"), i=CHROMSOME),
        fams = expand(os.path.join(SNPDB, "chr{i}.fam"), i=CHROMSOME),
    output:
         os.path.join(SNPDB, "PLINK/inbred.bed"),
         os.path.join(SNPDB, "PLINK/inbred.bim"),
         os.path.join(SNPDB, "PLINK/inbred.fam"),
    params:
        prefix = SNPDB
    shell:
        "plink --bfile {params.prefix}/chr1 --merge-list {input.merge} --make-bed --out {params.prefix}/PLINK/inbred"

rule rel_distance_matrix:
    input:
        os.path.join(SNPDB, "PLINK/inbred.bed"),
        os.path.join(SNPDB, "PLINK/inbred.bim"),
        os.path.join(SNPDB, "PLINK/inbred.fam"),
    output:
        os.path.join(SNPDB,"PLINK/inbred_grm.rel"),
        os.path.join(SNPDB,"PLINK/inbred_grm.rel.id"),
    params:
        prefix = os.path.join(SNPDB, "PLINK")
    threads: 12
    shell:
        "plink --bfile {params.prefix}/inbred " # need .bed, .bim, .fam
        "--make-rel square "
        "--out {params.prefix}/inbred_grm "
        "--threads {threads}"

rule rel_dist_ghmap:
    input:
        rel = os.path.join(SNPDB,"PLINK/inbred_grm.rel"),
        names = os.path.join(SNPDB,"PLINK/inbred_grm.rel.id"),
    output:
        rel = GENETIC_REL,
    run:
        samples = []
        with open(input.names, 'r') as name:
            for n in name:
                samples.append(n.strip().split("\t")[-1])
        samples = "\t".join(samples)
        with open(input.rel, 'r') as dist:
            matrix = dist.read()
        with open(output.rel, 'w') as out:
            out.write("#"+samples+"\n")
            out.write(matrix)

######################### Prepare Variant Annotation File #########################################
rule unGZip:
    input: os.path.join(VEP_DIR, "chr{i}.pass.vep.txt.gz"),
    output: temp(os.path.join(VEP_DIR, "chr{i}.pass.vep.txt"))
    shell:
        "zcat {input} > {output}"
        
rule annotateSNPs:
    input: 
        vep = os.path.join(VEP_DIR, "chr{i}.pass.vep.txt"),
    output: os.path.join(SNPDB, "chr{i}.snp.annot.txt")
    params:
        bin = HBCGM_BIN,
    shell:
        "{params.bin}/haplomap annotate -t snp -o {output} {input.vep} "


######################## START TO RUN Main Program #####################################
# find haplotypes
rule eblocks:
    input: 
        snps = os.path.join(SNPDB, "chr{i}.snp.txt"),
        gene_anno =  os.path.join(SNPDB, "chr{i}.snp.annot.txt"),
        strains = "MPD_{ids}/trait.{ids}.txt",
    output: 
        hb = protected("MPD_{ids}/chr{i}.snp.hblocks.txt"),
        snphb = protected("MPD_{ids}/chr{i}.snp.haplotypes.txt")
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
        hb = "MPD_{ids}/chr{i}.snp.hblocks.txt",
        trait = "MPD_{ids}/trait.{ids}.txt",
        gene_exprs = GENE_EXPRS,
        rel = GENETIC_REL
    output: "MPD_{ids}/chr{i}.snp.results.txt"
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

rule ghmap_aggregate:
    input: 
        res = ["MPD_{ids}/chr%s.snp.results.txt"%c for c in CHROMSOME]
    output: temp("MPD_{ids}_snp.results.txt")
    run:
        # read input
        dfs = []
        for p in input.res:
            case = pd.read_table(p, skiprows=5, dtype=str)
            dfs.append(case)
        result = pd.concat(dfs)
        # read header, first 5 row
        headers = []
        with open(p, 'r') as r:
            for i, line in enumerate(r):
                headers.append(line)
                if i >= 4: break 

        if os.path.exists(output[0]): os.remove(output[0])
        # write output
        with open(output[0], 'a') as out:
            for line in headers:
                out.write(line)
            ## Table
            result.to_csv(out, sep="\t", index=False)

rule mesh:
    input: 
        res = expand("MPD_{ids}_snp.results.txt", ids=IDS),
        json = MPD2MeSH,
    output: expand("MPD_{ids}_snp.results.mesh.txt", ids=IDS)
    params: 
        # mesh = lambda wildcards: MESH_DICT[wildcards.ids]
        gnnhap = config['HBCGM']['GNNHAP'],
        res_dir = config['HBCGM']['WORKSPACE'],
        bundle = config['HBCGM']['GNNHAP_BUNDLE'],
    threads: 24
    shell:
        "python {params.gnnhap}/GNNHap/predict.py "
        "--bundle {params.bundle} " 
        "--hbcgm_result_dir {params.res_dir} "
        "--mesh_terms {input.json} --num_cpus {threads} "


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


