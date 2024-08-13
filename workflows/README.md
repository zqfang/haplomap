# Snakemake pipelines for Variant calling

## Usage
### 1. Edit the `config.yaml` file for required files:
Install `snakemake`, `ensembl-vep` first.

update the path to:
- BAM_DIR
- GENOME
- VEP

And the output directory:
- BCFTOOL (bcftools pipeline)
  - WORKSPACE

- GATK (gatk pipeline)
  - WORKSPACE

### 2 Run on local computer
```shell
# modify the file path in haplomap and run with 12 cores
snakemake -s bcftools.call.smk  --configfile config.yaml \
          -k -p -j 12   
```

or 
```shell
# modify the file path in haplomap and run with 12 cores
snakemake -s gatk.call.smk  --configfile config.yaml \
          -k -p -j 12   
```

## About variant calling pipeline:
1. bcftools call 
  - prefered pipeline for inbred mouse and haplomap input

2. GATK best practice
  - designed for human genetics 
  - users are responsible for tuning ``VQSR`` or ``hardfiltering`` parameters if use non-human data


**Caution !**: Both pipelines take a long time to run.


## Ensemble VEP

### Install VEP

[Download and Install](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html)

or use conda to install `vep`

### Install Genome
Install data for offline mode
```shell
INSTALL.pl -a cfp -s mus_musculus -y GRCm38 --CONVERT --PLUGINS CADD,GO,TSSDistance,LoF,SpliceAI
```

or Human
```shell
INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 --CONVERT
```


### Annotation
e.g
```shell
bcftools view -f .,PASS ${input.vcf} | \
        ${params.VEPBIN}/vep --fasta ${input.reference} ${params.genome_build} \
        --format vcf --fork ${threads} --hgvs --force_overwrite \
        --uniprot --domains --symbol --regulatory --distance 1000 --biotype \
        --gene_phenotype MGI --check_existing  --pubmed --numbers \
        --offline --cache --variant_class \
        --gencode_basic --no_intergenic --individual all \
        -o ${output} --tab --compress_output gzip \
```

### Notes


## Why not GATK for inbred population ?

One of my colleague who studies mouse genetics, said, 

> I tried the haplotype caller from GATK. But it seems that the haplotype caller is designed for heterogeneous genome like human than for mice. Therefore, the result coming out of HC is worse than samtools, as I manually inspected a few regions that HC calls didn't make sense.

> In addition, in one of their mouse genomic paper that we reviewed, they even skipped the second recalibration step. We asked them why and they said it was because of the same reason: good for human but not that good for the homogeneous inbred mouse.


However, my experiences with GATK (>v4.0) is as good as bcftools.