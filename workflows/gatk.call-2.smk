import os
from snakemake.shell import shell
############### Globals ########################
WORKSPACE = "/data/bases/fangzq"
workdir: WORKSPACE

GENOME="/home/fangzq/genome/mouse/GRCm38_68.fa"
dbSNP="/home/fangzq/genome/mouse/mgp.v5.merged.snps_all.dbSNP142.sorted.vcf"
dbINDEL="/home/fangzq/genome/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf"
STRAINS_FILE = "/data/bases/fangzq/strain"
BAM_DIR = "/data/bases/fangzq/strains"
TMPDIR = "/home/fangzq/TMPDATA"

#CHROMSOME = [ str(c) for c in range(1,20)] + ["X", "Y", "MT"]
CHROMSOME = ['1'] + [ str(c) for c in range(10,20)] + [ str(c) for c in range(2,10)]+ ["MT", "X", "Y"]

# with open(STRAINS_FILE, 'r') as s:
#     STRAINS = s.read().strip().split()
 
STRAINS = ['129P2', '129S1', '129S5', 'AKR', 'A_J', 'B10', 
        'BPL', 'BPN', 'BTBR', 'BUB', 'B_C', 'C3H', 'C57BL10J',
        'C57BL6NJ', 'C57BRcd', 'C57LJ', 'C58', 'CBA', 'CEJ', 
        'DBA', 'DBA1J', 'FVB', 'ILNJ', 'KK', 'LGJ', 'LPJ', 
        'MAMy', 'MRL', 'NOD', 'NON', 'NOR', 'NUJ', 'NZB', 'NZO', 'NZW', 
        'PJ', 'PLJ', 'RFJ', 'RHJ', 'RIIIS', 'SEA', 'SJL', 'SMJ', 'ST', 'SWR', 'TALLYHO', 'RBF'] + \
        ['CAST', 'MOLF', 'PWD','PWK', 'SPRET', 'WSB']  # <- wild derived except MRL
# OUTPUT
VCF_VQSR = expand("VCFs/combined.chr{i}.VQSR.vcf.gz", i=CHROMSOME)
VCF_HFILTER = expand("VCFs/combined.chr{i}.hardfilter.vcf.gz", i=CHROMSOME)
VCF_HFILTER_PASS = expand("VCFs/combined.chr{i}.hardfilter.pass.vcf.gz", i=CHROMSOME)
VCF_RAW = expand("VCFs/combined.chr{i}.raw.vcf", i=CHROMSOME)
GVCF = expand("GVCF/{sample}.raw.g.vcf", sample=STRAINS)
# SNPs = expand("SNPs/combined.chr{i}.txt", i=CHROMSOME)


############## Rules ##########################
rule all:
    input: VCF_HFILTER_PASS, #VCF_VQSR

# include: "rules/gatk.getbam.smk"

rule sampleCalling:
    input:
        dbSNP = dbSNP,
        genome = GENOME,
        genomedict = GENOME.replace(".fa", ".dict"),
        #bam = "BAM/{sample}.marked.fixed.BQSR.bam",
        bam = os.path.join(BAM_DIR, "{sample}/output.GATKrealigned.Recal.bam")
    output: 
        # save for future use when more strains are added
        gvcf= protected("GVCF/{sample}.raw.g.vcf"), 
        gvcfi= protected("GVCF/{sample}.raw.g.vcf.idx") 
    threads: 4
    log: "logs/{sample}.haplotypecaller.log"
    params:
        java_ops="-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    shell:
        "gatk --java-options '{params.java_ops}' HaplotypeCaller "
        "-ERC GVCF "
        "--native-pair-hmm-threads {threads} "
        "--dbsnp {input.dbSNP} "
        "-R {input.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} 2> {log} "

## temp output for combineGVCFs.
rule chroms: 
    output: temp(expand("GVCF/chr{i}.combine.tmp", i = CHROMSOME))
    run:
        for out in output:
            shell("touch %s"%out)

# CombinedGVCFs, has to be an interval 
## parallel computing for joint_calling
# WARNING: This step used a lot of memory when combineGVCFs.        
rule combineGVCFs:
    input:
        genome=GENOME,
        chrom="GVCF/chr{i}.combine.tmp",
        gvcf=expand("GVCF/{sample}.raw.g.vcf", sample=STRAINS),
        gvcfi=expand("GVCF/{sample}.raw.g.vcf.idx", sample=STRAINS)
    output:
        gvcf=temp("GVCF/combined.chr{i}.g.vcf"),
        gvcfi=temp("GVCF/combined.chr{i}.g.vcf.idx")
    params:
        chrom=CHROMSOME,
        gvcf=" --variant ".join(expand("GVCF/{sample}.raw.g.vcf", sample=STRAINS)),
        java_ops="-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    log: "logs/chr{i}.combineGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' CombineGVCFs -L {wildcards.i} -R {input.genome} "
        "--variant {params.gvcf} -O {output.gvcf} 2> {log}"

# NOTE：GATK 相比于samtools 更容易找到杂合突变
rule jointCalling:
    input: 
        gvcf="GVCF/combined.chr{i}.g.vcf",
        gvcfi="GVCF/combined.chr{i}.g.vcf.idx",
        genome=GENOME
    output: 
        vcf="VCFs/combined.chr{i}.raw.vcf",
        vcfi="VCFs/combined.chr{i}.raw.vcf.idx",
    params:
        tmpdir=TMPDIR,
        java_ops= "-Xmx32G -Djava.io.tmpdir=%s"%TMPDIR
    log: "logs/chr{i}.GenotypeGVCFs.log"
    shell:
        "gatk --java-options '{params.java_ops}' "
        "GenotypeGVCFs -R {input.genome} -V {input.gvcf} -O {output.vcf} 2> {log}"


################# VQSR ###############################
# for human, first choice is VQSR
# for non-human, better to do is hardfilering
rule variantRecalSNPs:
    input: 
        vcf="VCFs/combined.chr{i}.raw.vcf",
        vcfi="VCFs/combined.chr{i}.raw.vcf.idx",
        genome=GENOME,
    output: 
        recal="VCFs/combined.chr{i}.snp.Recal",
        tranches= "VCFs/combined.chr{i}.snp.tranches",
        rscript= "VCFs/combined.chr{i}.snp.R", 
    params:
        dbsnp=f"dbsnp,known=true,training=true,truth=true,prior=2.0:{dbSNP}",  
    log: "logs/combined.chr{i}.snp.Recal.log"   
    shell:
        "gatk VariantRecalibrator -mode SNP "
        "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP "
        "--resource {params.dbsnp} "
        "-R {input.genome} -V {input.vcf} -O {output.recal} "
        "--tranches-file {output.tranches} --rscript-file {output.rscript} 2> {log}"

rule applyVQSR4SNPs:
    input:
        genome=GENOME,
        vcf="VCFs/combined.chr{i}.raw.vcf",
        vcfi="VCFs/combined.chr{i}.raw.vcf.idx",
        recal="VCFs/combined.chr{i}.snp.Recal",
        tranches= "VCFs/combined.chr{i}.snp.tranches",
        rscript= "VCFs/combined.chr{i}.snp.R",        
    output:
        vcf=temp("VCFs/combined.chr{i}.snp.VQSR.vcf"),
        vcfi=temp("VCFs/combined.chr{i}.snp.VQSR.vcf.idx"),
    log: "logs/combined.chr{i}.snp.VQSR.log"
    shell:
        "gatk ApplyVQSR -mode SNP "
        "-R {input.genome} -V {input.vcf} -O {output} "
        "--truth-sensitivity-filter-level 99.5 "
        "--tranches-file {input.tranches} "
        "--rscript-file {input.rscript} "
        "--recal-file {input.recal} 2> {log} "

rule variantRecalINDELs:
    input: 
        vcf="VCFs/combined.chr{i}.raw.vcf",
        vcfi="VCFs/combined.chr{i}.raw.vcf.idx",
        genome=GENOME,
    output: 
        recal ="VCFs/combined.chr{i}.indel.Recal",
        tranches = "VCFs/combined.chr{i}.indel.tranches",
        rscript = "VCFs/combined.chr{i}.indel.plots.R",
    params:
        #format: "{},known={},training={},truth={},prior={}:{}"
        dbindel=f"snps_all,known=true,training=true,truth=true,prior=12.0:{dbINDEL}",
        dbsnp=f"indels,known=true,training=true,truth=true,prior=2.0:{dbSNP}", 
    log: "logs/combined.chr{i}.indel.Recal.log"
    shell:
        "gatk VariantRecalibrator -mode INDEL "
        "--max-gaussians 4"
        "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum "
        "--resource {params.dbindel} "
        "--resource {params.dbsnp} "
        "-R {input.genome} -V {input.vcf} -O {output.recal} "
        "--tranches-file {output.tranches} --rscript-file {output.rscript} 2> {log}"

rule applyVQSR4INDEL:
    input:
        genome = GENOME,
        vcf = "VCFs/combined.chr{i}.raw.vcf",
        vcfi= "VCFs/combined.chr{i}.raw.vcf.idx",
        recal = "VCFs/combined.chr{i}.indel.Recal",
        tranches = "VCFs/combined.chr{i}.indel.tranches",
        rscript = "VCFs/combined.chr{i}.indel.plots.R", 
    output:
        vcf=temp("VCFs/combined.chr{i}.indel.VQSR.vcf"),
        vcfi=temp("VCFs/combined.chr{i}.indel.VQSR.vcf.idx"),
    log: "logs/combined.chr{i}.indel.VQSR.log"
    shell:
        "gatk ApplyVQSR -mode INDEL "
        "-R {input.genome} -V {input.vcf} -O {output.vcf} "
        "--tranches-file {input.tranches} "
        "--rscript-file {input.rscript} "
        "--recal-file {input.recal} "
        "--truth-sensitivity-filter-level 99.0 2> {log} "

rule mergeVQSRVCFs:
    input:
        snp = "VCFs/combined.{chrom}.snp.VQSR.vcf",
        snpi = "VCFs/combined.{chrom}.snp.VQSR.vcf.idx",
        indel = "VCFs/combined.{chrom}.indel.VQSR.vcf",
        indeli = "VCFs/combined.{chrom}.indel.VQSR.vcf.idx",
    output:
        vcf=protected("VCFs/combined.{chrom}.VQSR.vcf"),
        vcfi=protected("VCFs/combined.{chrom}.VQSR.vcf.idx")
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} "
        "-O {output.vcf} 2>/dev/null"

rule compressVQSRVCFs:
    input: "VCFs/combined.{chrom}.VQSR.vcf"
    output: protected("VCFs/combined.{chrom}.VQSR.vcf.gz")
    shell:
        """bgzip -f {input} 
           tabix -p vcf {output}
        """

################# Hard filering ######################
rule selectSNPs:
    input: 
        vcf="VCFs/combined.{chrom}.raw.vcf",
        vcfi="VCFs/combined.{chrom}.raw.vcf.idx",
    output: 
        vcf=temp("VCFs/combined.{chrom}.snp.vcf"),
        vcfi=temp("VCFs/combined.{chrom}.snp.vcf.idx"),
    shell:
        "gatk SelectVariants -select-type SNP " 
        "-V {input.vcf} -O {output.vcf} 2>/dev/null"

rule hardFilterSNPs:
    input: 
        vcf="VCFs/combined.{chrom}.snp.vcf",
        vcfi="VCFs/combined.{chrom}.snp.vcf.idx"
    output: 
        vcf="VCFs/combined.{chrom}.snp.filter.vcf",
        vcfi="VCFs/combined.{chrom}.snp.filter.vcf.idx",
    log: "logs/combined.{chrom}.snp.filter.log"
    shell:
        "gatk VariantFiltration "
        "-filter 'QD < 2.0' --filter-name 'QD2' " 
        "-filter 'QUAL < 30.0' --filter-name 'QUAL30' " 
        "-filter 'SOR > 3.0' --filter-name 'SOR3' " 
        "-filter 'FS > 60.0' --filter-name 'FS60' " 
        "-filter 'MQ < 40.0' --filter-name 'MQ40' " 
        "-filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' " 
        "-filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' " 
        "-V {input.vcf} -O {output.vcf} 2> {log}"    

rule selectINDELs:
    input: 
        vcf="VCFs/combined.{chrom}.raw.vcf",
        vcfi="VCFs/combined.{chrom}.raw.vcf.idx",
    output: 
        vcf=temp("VCFs/combined.{chrom}.indel.vcf"),
        vcfi=temp("VCFs/combined.{chrom}.indel.vcf.idx"),
    shell:
        "gatk SelectVariants -select-type INDEL " 
        "-V {input.vcf} -O {output.vcf} 2>/dev/null"   

rule hardFilterINDELs:
    input: 
        vcf="VCFs/combined.{chrom}.indel.vcf",
        vcfi="VCFs/combined.{chrom}.indel.vcf.idx",
    output: 
        vcf=temp("VCFs/combined.{chrom}.indel.filter.vcf"),
        vcfi=temp("VCFs/combined.{chrom}.indel.filter.vcf.idx"),
    log: "logs/combined.{chrom}.indel.filter.log"
    shell:
        "gatk VariantFiltration " 
        "-filter 'QD < 2.0' --filter-name 'QD2' " 
        "-filter 'QUAL < 30.0' --filter-name 'QUAL30' " 
        "-filter 'FS > 200.0' --filter-name 'FS200' " 
        "-filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' " 
        "-V {input.vcf} -O {output.vcf} 2> {log}"     

rule mergeHardVCFs:
    input:
        snp= "VCFs/combined.{chrom}.snp.filter.vcf",
        snpi= "VCFs/combined.{chrom}.snp.filter.vcf.idx",
        indel= "VCFs/combined.{chrom}.indel.filter.vcf",
        indeli= "VCFs/combined.{chrom}.indel.filter.vcf.idx",
    output: 
        vcf=protected("VCFs/combined.{chrom}.hardfilter.vcf"),
        vcfi=protected("VCFs/combined.{chrom}.hardfilter.vcf.idx")
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} "
        "-O {output.vcf} 2>/dev/null"

rule compressHardVCF:
    input: "VCFs/combined.{chrom}.hardfilter.vcf"
    output: protected("VCFs/combined.{chrom}.hardfilter.vcf.gz")
    shell:
        """bgzip -f {input} 
           tabix -p vcf {output}
        """

rule extractPASS:
    input: 
        vcf="VCFs/combined.{chrom}.hardfilter.vcf",
        genome=GENOME
    output: "VCFs/combined.{chrom}.hardfilter.pass.vcf.gz"
    shell:
        "gatk SelectVariants -R {input.genome} -V {input.vcf} -O {output} "
        " -select 'vc.isNotFiltered()' 2>/dev/null"