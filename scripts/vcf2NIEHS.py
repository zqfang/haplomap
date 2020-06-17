import time, sys, os, glob, gzip
from collections.abc import Iterable

def vcf2niehs(invcf, outdir, chromosome, qual_samtools=50, heterzygote_cutoff = 20):
    
    ## strains
    # if isinstance(strains, str) or not isinstance(strains, Iterable):
    #     with open(strains, 'r') as ss:
    #         strains = ss.read().strip().split()
    # # rename B_C to BALB
    # if strains.index("B_C") != -1:
    #     strains[strains.index("B_C")] = "BALB"
    os.makedirs(outdir, exist_ok=True)
    ##### filtering parameters
    qualCutoffForSamtools = qual_samtools
    cutoffForPLscoreForNonHomozygousAlt = heterzygote_cutoff

    numOfAlt, numGoodAlt, theGoodAlt = 0, 0, 0
    totalVariant, totalInterested, numNonPass, numINDEL, numLowQual, numMultiAlt, numGoodSNP = 0, 0, 0, 0, 0, 0, 0
    nucleotide = {"A": 1, "C":1, "G":1, "T":1} 
    
    # prepared output
    if isinstance(chromosome, str) or not isinstance(chromosome, Iterable):
        chromosome = [chromosome.strip("chr")]
    outputdict = { str(k) : open(os.path.join(outdir, f"chr{k}.full.txt" ), 'w') for k in chromosome }
    output_compact = { str(k) : open(os.path.join(outdir, f"chr{k}.txt" ), 'w') for k in chromosome }
    # add header
    # for chrom, output in outputdict.items():
    #     output.write("LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES\t")
    #     output.write("\t".join(strains)+ "\n")
    # for chrom, output in output_compact.items():
    #     output.write("\t".join(["C57BL/6J"] + strains)+ "\n")
    # parse VCF
    # CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT	SAMPLES1 ...
    lineCount = 0
    sname = os.path.basename(invcf)
    
    if invcf.endswith("gz"):
        vcf = gzip.open(invcf, 'rt') 
    else:
        vcf = open(invcf, 'r')
    strains = []
    for line in vcf:
        # skip description lines
        if line.startswith("##"): continue
        # write header 
        if line.startswith("#CHROM"): 
            strain_name = line.split("FORMAT")[-1]
            strains += strain_name.strip().split("\t")
            for chrom, output in outputdict.items():
                output.write("LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES\t")
                output.write("\tC57BL/6J"+strain_name)
            for chrom, output in output_compact.items():
                output.write("\tC57BL/6J"+strain_name)
            continue

        # parse data
        newline = line.strip().split("\t")

        lineCount +=1
        totalVariant +=1            
        ## check vcf format
        if len(newline) < 9: 
            sys.exit("Error! format not recognized")
        # skip chr if not interested 
        if newline[0] not in chromosome: continue
        totalInterested +=1
        
        # FIXME: gatk not pass strings: MQ40, SOR3 ... ?
        # gatk -> [PASS, LowQual, '.' ], 
        # bcftools -> [PASS, LowQual, '.' ],
        # -> not filtering has been done when it's '.'
        # if newline[6] in ["PASS", ".", "INDEL"]:
        #     numNonPass +=1
        #     continue
        ### first, get the ref/alternative alleles

        # NOTE: To remove, GATK INDELs
        # run gatk SelectVariant --select-type SNP...
        ref = newline[3]
        # skip indels
        # -1 if could not find pattern "INDEL"
        if len(ref) > 1 or newline[7].find("INDEL") != -1: 
            numINDEL +=1
            continue
        
        # QUAL: The Phred-scaled probability that a REF/ALT polymorphism exists.
        # QUAL: -10 * log (1-p)
        # it's not a very useful property for evaluating the quality of a variant call 
        # FIXME: the cutoff for GATK (30 default) -> VQSR or Hardfilering first!
        if float(newline[5]) < qualCutoffForSamtools:
            numLowQual +=1
            continue

        ###find the entries for GT & PL
        # GATK
        # if newline[8] == 'GT:AD:DP:GQ:PL':
        #     GTind, PLind = 0, 4
        # elif newline[8] == 'GT:AD:DP:GQ:PGT:PID:PL':
        #     GTind, PLind = 0, 6
        # # samtools
        # elif newline[8] == "GT:PL:DP:DV:SP:DP4:DPR:GP:GQ":
        #     GTind, PLind = 0, 1
        # else:
        IDS = newline[8].split(":") 
        GTind, PLind = -1, -1
        try:
            GTind = IDS.index('GT')
            PLind = IDS.index('PL')
        except ValueError:
            sys.exit("Error! PL or GT NOT Found!") 
        
        chrom = newline[0]
        pos = newline[1]  
        alts = newline[4].split(",")
        numOfAlt = len(alts)
        # FIXME: 
        # if numOfAlt > 1: sys.exit("Bad Alt")

        hasAlt= [0] * len(alts)
        alleles = [-1] * len(strains)
        # CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT	AKR
        for s, strain in enumerate(strains):
            ## identify the allele for each strai
            strainFormats = newline[9+s].split(":")

            # if "" -> NN
            if strainFormats[GTind] == "./." or strainFormats[GTind] == ".|." or (not strainFormats[GTind]) or (not strainFormats[PLind]):         
                alleles[s] = "NN"
                continue

            GTs = strainFormats[GTind].replace("|", "/").split("/")
            if len(GTs) != 2:
                sys.exit(f"Improper GT: {strainFormats[GTind]} in {line}\n")

            if GTs[0] == "." or GTs[0] != GTs[1]:
                alleles[s] = "NN" # nocall/heterozyte -> NN
                continue
            ## now this is supposed to be a homogyzous call. Check whether the call is of good quality
            # find GT       
            GTs = [int(s) for s in GTs]
            index = (GTs[0]+ 1) * (GTs[0] + 2) // 2 - 1 # 0,2
            PLs = strainFormats[PLind].split(",")
            # PLs -> ['r/r', 'r/a', 'a/a']
            PLs = [float(s) for s in PLs]
            if len(PLs) != ((numOfAlt + 1) * (numOfAlt + 2) // 2):
                sys.exit(f"Improper PL found: {strainFormats[PLind]}\n{line}")
            
            minScore = 10000000000
            for i, pl in enumerate(PLs):
                if i == index: continue 
                #  PL_het - PL_homo 
                if (PLs[i] - PLs[index]) < minScore:
                    minScore = PLs[i] - PLs[index] 
            if minScore >= cutoffForPLscoreForNonHomozygousAlt:
                #  PL_het - PL_homo >= 20 
                # log10 (pHe / pHo) <= -2
                # pHe / pHo <= 0.01
                ## this is a good call
                if GTs[0] > 0:
                    # select alt
                    alleles[s] = f"{alts[GTs[0]-1]}{alts[GTs[0]-1]}"
                    hasAlt[GTs[0]-1] = 1
                else:
                    alleles[s] = f"{ref}{ref}"
            else:
                alleles[s] = "NN"

        numGoodAlt = 0
        theGoodAlt = -1
        for i, k in enumerate(hasAlt):
            if k == 1:
                numGoodAlt +=1
                theGoodAlt = i
        if numGoodAlt == 1:
            if theGoodAlt == -1:
                sys.exit("Interal error: the good Alt not found properly")
            if alts[theGoodAlt] not in nucleotide.keys():
                numINDEL +=1

            ## this is a good SNP across the strains. save it for output 
            allAlleles = "\t".join(alleles)
            alt = alts[theGoodAlt]
            # write output here
            ## header 
            # ## LOCAL_IDENTIFIER SS_ID CHROMOSOME ACCESSION_NUM POSITION STRAND ALLELES + strains ...
            D = f"SNP_{chrom}_{pos}"
            outline = f"{D}\t{D}\t{chrom}\t{D}\t{pos}\t.\t{ref}/{alt}\t{allAlleles}\n"
            outputdict[chrom].write(outline)
            alleles_pattern = [ref] + [s[0] for s in alleles]
            alleles_compact = ''.join(alleles_pattern).replace("N", "?")
            outline_compact = f"{D}\t{chrom}\t{pos}\t{alleles_compact}\n"
            output_compact[chrom].write(outline_compact)
            numGoodSNP += 1
        else:
            if numGoodAlt == 0:
                numLowQual += 1
            else:
            ## multiple good alternatives found;
                numMultiAlt +=1  
        
        if lineCount % 100000 == 0:
            print(f"Sample: {sname} - Finished line: {lineCount}", file=sys.stderr)        

    # close file
    for chrom, output in outputdict.items(): output.close()
    for chrom, output in output_compact.items(): output.close()
    vcf.close()

    print(f"Sample: {sname} - total: {totalVariant}, interested: {totalInterested}, notPass: {numNonPass}, " +\
          f"indels: {numINDEL}, lowQual: {numLowQual}, multiAlt: {numMultiAlt}, good: {numGoodSNP}", file=sys.stderr)


vcf2niehs(snakemake.input['vcf'], snakemake.params['outdir'], 
          snakemake.params['chrom'],
          snakemake.params['qual_samtools'], 
          snakemake.params['heterzygote_cutoff'])


# if __name__ == '__main__':
#     # inputs and outputs
#     chromosome = list(range(1, 20)) + [ "X", "Y"]
#     for i in chromosome:
#         invcf = f"VCFs/combined.chr{i}.raw.vcf.gz"
#         outdir = f"SNPs"
#         vcf2niehs(invcf, outdir, str(i), 50,  20)