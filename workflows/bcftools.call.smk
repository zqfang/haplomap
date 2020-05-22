
import os
from snakemake.shell import shell
############### Globals ########################
WORKSPACE = "/data/bases/fangzq/20200429"
workdir: WORKSPACE

GENOME="/home/fangzq/genome/mouse/GRCm38_68.fa"
dbSNP="/home/fangzq/genome/mouse/mgp.v5.merged.snps_all.dbSNP142.sorted.vcf"
dbINDEL="/home/fangzq/genome/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf"
STRAINS_FILE = "/data/bases/fangzq/strain"
BAM_DIR = "/data/bases/fangzq/strains"

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
#CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["X", "Y"]
CHROMSOME = ['4','Y']
# with open(STRAINS_FILE, 'r') as s:
#     STRAINS = s.read().strip().split()
 
STRAINS = ['129P2', '129S1', '129S5', 'AKR', 'A_J', 'B10', 
        'BPL', 'BPN', 'BTBR', 'BUB', 'B_C', 'C3H', 'C57BL10J',
        'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
        'DBA', 'DBA1J', 'FVB', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
        'MAMy', 'MRL','NOD', 'NON', 'NOR', 'NUJ', 'NZB', 'NZO', 'NZW', 
        'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR', 'TALLYHO', 'RBF'] + \
         ['CAST', 'MOLF', 'PWD','PWK', 'SPRET', 'WSB']  # <- wild derived except MRL
# OUTPUT

VCF_HFILTER = expand("VCFs/combined.chr{i}.hardfilter.vcf.gz", i=CHROMSOME)
VCF_HFILTER_PASS = expand("VCFs/combined.chr{i}.hardfilter.pass.vcf.gz", i=CHROMSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf", i=CHROMSOME)
VCF_STATS = expand("VCFs/combined.chr{i}.raw.vcf.gz.stats", i=CHROMSOME)

SNPDB = expand("SNPs/chr{i}.txt", i=CHROMSOME)

###########################################################################################
rule target:
    input: VCF_STATS, VCF_HFILTER_PASS,SNPDB

# samtools-bcftools-calling
rule faidx: 
    input: GENOME
    output: GENOME + ".fai"
    shell:
        "samtools faidx input.fa"

rule getChromSize:
    input: GENOME + ".fai"
    output: temp(expand("chr{i}.tmp",i=CHROMSOME))
    params:
        cmd='{print $1":1-"$2}',
        chrom=CHROMSOME
    #shell:
        # snakemake use bash -c mode 
        # Escape certain characters, such as \t by \\t, $ by \$, and { by {{.
        # Use triple quotation marks to surround the command line call.
        # """awk -F"\\t" '{{print $1":1-"$2}}' {input} | while read line;
        #    do echo -e $line > region.$line done
        # """
        # this works too
        # "awk -F'\\t' '{params.cmd}' {input} | while read line; "
        # "do echo -e $line > region.$line done"
    run:
        with open(input[0]) as cords:
            for c in cords:
                c1 = c.strip().split("\t")
                cmd = f"echo -e '{c1[0]}:1-{c1[1]}' > chr{c1[0]}.tmp"
                shell(cmd)
            
rule bcftools_call:
    input:
        genome = GENOME,
        bam = expand(os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam"), sample=STRAINS),
        #bam = expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS),
        chroms_region="chr{i}.tmp"
    output: 'VCFs/combined.chr{i}.raw.vcf'
    params:
        #bam = " ".join(expand("BAM/{sample}.marked.fixed.BQSR.bam", sample=STRAINS))
        bam=" ".join(expand(os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam"), sample=STRAINS))
    run:
        with open(input.chroms_region) as reg:
            region = reg.read().strip()
        cmd = "bcftools mpileup "+\
              "-a DP,AD,ADF,ADR,SP,INFO/AD "+\ 
              "-E -F0.25 -Q0 -p -m3 -d500 -r %s "%region +\
              "-Ou -f {input.genome} {params.bam} | " +\  
              "bcftools call -mv -f GQ,GP -Ov  > {output}"
        shell(cmd)

rule tabix:
    input: "VCFs/combined.{chr}.raw.vcf"
    output: 
        temp("VCFs/combined.{chr}.raw.vcf.gz"),
        temp("VCFs/combined.{chr}.raw.vcf.gz.tbi")
    shell:
        """bgzip -f {input} 
           tabix -p vcf {output[0]}
        """

rule bcftools_stats:
    input: 
        genome=GENOME,
        vcf="VCFs/combined.{chr}.raw.vcf.gz",
        vcfi="VCFs/combined.{chr}.raw.vcf.gz.tbi"
    output: 
        "VCFs/combined.{chr}.raw.vcf.gz.stats"
    shell:
        "bcftools stats -F {input.genome} -s - {input.vcf} > {output}"

rule bcftools_plot:
    input: "VCFs/combined.{chr}.raw.vcf.gz.stats"
    output: "figures/combined.{chr}.raw.vcf.gz.stats.pdf"
    shell:
        "plot-vcfstats -p figures -T {wildcards.chr} {input}"

rule bcfcall_filtering:
    input: 
        vcf="VCFs/combined.{chr}.raw.vcf.gz",
        vcfi="VCFs/combined.{chr}.raw.vcf.gz.tbi"
    output: 
        "VCFs/combined.{chr}.hardfilter.pass.vcf.gz"
    shell: 
        "bcftools filter -Oz -o {output} -s LOWQUAL -i'%QUAL>20' {input}"

rule unbigzip:
    input: "VCFs/combined.{chr}.hardfilter.pass.vcf.gz"
    output: temp("VCFs/combined.{chr}.hardfilter.pass.vcf")
    shell:
        "bgzip -d {input}"

rule vcf2strains:
    input:  
        "VCFs/combined.{chr}.hardfilter.pass.vcf"
    output: 
        temp("SNPs/{chr}.strains.temp")
    shell:
        # NOTE: '\t' is default delim for cut
        "head -n 1000 {input} | grep '^#CHROM' | "
        "cut -f10-  > {output}"  
    
rule vcf2niehs:
    input:  
        # vcf = "VCFs/combined.chr{i}.raw.vcf", 
        vcf = "VCFs/combined.chr{i}.hardfilter.pass.vcf",
        strains = "SNPs/chr{i}.strains.temp"
    output: 
        protected("SNPs/chr{i}.txt")
    params:
        outdir= "SNPs",
        chrom="{i}",
        qual_samtools=50, 
        heterzygote_cutoff = 20
    script:
        "../scripts/vcf2NIEHS.py"