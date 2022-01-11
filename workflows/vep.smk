



import os
############### Globals ########################
# configfile: "config.yaml"
workdir: config['BCFTOOLS']['WORKSPACE']

GENOME = config['GENOME']
dbSNP = config['dbSNP']
BAM_DIR = config['BAM_DIR']
STRAINS = sorted(config['STRAINS'])

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X"]

# OUTPUT
VCF_PASS = expand("VCFs/chr{i}.vcf.gz", i=CHROMSOME)
VEP_INDEL = expand("VEP/chr{i}.indel.pass.vep.txt.gz", i=CHROMSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf.gz", i=CHROMSOME)
VCF_STATS = expand("VCFs/combined.chr{i}.raw.vcf.gz.stats", i=CHROMSOME)
ANNO_EBLOCK = expand("INDELs_Annotation/chr{i}.annotation.eblocks.txt", i=CHROMSOME)

###########################################################################################
rule target:
    input: ANNO_EBLOCK, #VCF_PASS, VEP_INDEL, ANNO_EBLOCK # SNPDB,VCF_RAW, VCF_PASS, SNPDB, #VEP_ANNO, #VCF_STATS



# rule VEP4Indels:
#     """emsemble-vep"""
#     input: 
#         vcf="VCFs/{chr}.indel.vcf.gz",
#         reference=GENOME,
#     output: "VEP/{chr}.indel.pass.vep.txt.gz"
#     params:
#         #genome_build = " -a GRCm38 --species mus_musculus ",
#         genome_build = config['VEP']['GENOME_BUILD'],
#         VEPBIN = config['VEP']['BIN'],
#         extra=" --dir_cache "  + config['VEP']['CACHE_DIR']
#     threads: 4
#     shell:
#         ## emsemble-vep
#         # https://github.com/Ensembl/ensembl-vep
#         "bcftools view -f .,PASS {input.vcf} -v indels | "
#         "{params.VEPBIN}/vep --fasta {input.reference} {params.genome_build} "
#         "--format vcf --fork {threads} --hgvs --force_overwrite "
#         "--uniprot --domains --symbol --regulatory --distance 1000 --biotype "
#         "--gene_phenotype MGI --check_existing  --pubmed --numbers "
#         "--offline --cache --variant_class "
#         "--gencode_basic --no_intergenic --individual all "
#         "-o {output} --tab --compress_output gzip"

rule Indel2anno:
    input:
        vcf="VCFs/{chr}.indel.vcf.gz",
        vep="VEP/{chr}.indel.pass.vep.txt.gz",
    output:
        anno="INDELs_Annotation/{chr}.annotation.txt",
        eblock="INDELs_Annotation/{chr}.annotation.eblocks.txt",
    shell:
        "python /home/fangzq/github/HBCGM/scripts/annotateINDEL.py {input.vcf} {input.vep} {output.anno} {output.eblock}"