import time, sys, os, glob, gzip
from collections.abc import Iterable

class VCF:
    def __init__(self, invcf: str, output: str, qual=50, heterozyote_thresh=20, samples=None):
        self.vcf = invcf
        self.output = output
        self.qual = qual
        self.het = heterozyote_thresh
        self.sname = os.path.basename(invcf)
        self.nucleotide = {"A": 1, "C":1, "G":1, "T":1} 
        # os.makedirs(self.outdir, exist_ok=True)
        self.samples = samples
        self.outputdict = open(output.replace(".txt", ".full.txt" ), 'w') 
        self.output_compact = open(output, 'w')         
        if invcf.endswith("gz"):
            self._vcf = gzip.open(invcf, 'rt') 
        else:
            self._vcf = open(invcf, 'r')
        
    def __del__(self):
        # close file
        self.outputdict.close()
        self.output_compact.close()
        self._vcf.close()

    def toNIEHS(self):
        ##### filtering parameters
        numOfAlt, numGoodAlt, theGoodAlt = 0, 0, 0
        totalVariant, totalInterested, numNonPass, numINDEL, numLowQual, numMultiAlt, numGoodSNP = 0, 0, 0, 0, 0, 0, 0
        lineCount = 0
        for line in self._vcf:
            # skip description lines
            if line.startswith("##"): continue
            # write header 
            if line.startswith("#CHROM"): 
                strain_name = line.split("FORMAT")[-1]
                strain_name = strain_name.replace("B_C", "BALB")
                #strain_name = strain_name.replace("A_J", "A/J")
                self.strains = strain_name.strip().split("\t")
                self.outputdict.write("LOCAL_IDENTIFIER\tSS_ID\tCHROMOSOME\tACCESSION_NUM\tPOSITION\tSTRAND\tALLELES\t")
                if self.samples is not None:
                    strain_name = "\t" + "\t".join(self.samples) + "\n"
                self.outputdict.write("C57BL/6J"+strain_name)
                self.output_compact.write("C57BL/6J"+strain_name)
                continue
            # parse data
            newline = line.strip().split("\t")

            lineCount +=1
            totalVariant +=1            
            ## check vcf format
            if len(newline) < 9: 
                sys.exit("Error! format not recognized")
            # skip chr if not interested 
            # if newline[0] not in chromosome: continue
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
            if float(newline[5]) < self.qual:
                numLowQual +=1
                continue

            ###find the entries for GT & PL
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
            alleles = [-1] * len(self.strains)
            # CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT
            for s, strain in enumerate(self.strains):
                ## identify the allele for each strain
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
                if minScore >= self.het:
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
                if alts[theGoodAlt] not in self.nucleotide.keys():
                    numINDEL +=1

                    
                ## this is a good SNP across the strains. save it for output 
                if self.samples is not None:
                    alleles = [alleles[self.strains.index(s)] for s in self.samples]

                allAlleles = "\t".join(alleles)
                alt = alts[theGoodAlt]
                # write output here
                ## header 
                # ## LOCAL_IDENTIFIER SS_ID CHROMOSOME ACCESSION_NUM POSITION STRAND ALLELES + strains ...
                D = f"SNP_{chrom}_{pos}"
                outline = f"{D}\t{D}\t{chrom}\t{D}\t{pos}\t.\t{ref}/{alt}\t{allAlleles}\n"
                self.outputdict.write(outline)
                alleles_pattern = [ref] + [s[0] for s in alleles]
                alleles_compact = ''.join(alleles_pattern).replace("N", "?")
                outline_compact = f"{D}\t{chrom}\t{pos}\t{alleles_compact}\n"
                self.output_compact.write(outline_compact)
                numGoodSNP += 1
            else:
                if numGoodAlt == 0:
                    numLowQual += 1
                else:
                ## multiple good alternatives found;
                    numMultiAlt +=1  
            
            if lineCount % 100000 == 0:
                print(f"Sample: {self.sname} - Finished line: {lineCount}", file=sys.stderr)        

        print(f"Sample: {self.sname} - total: {totalVariant}, interested: {totalInterested}, notPass: {numNonPass}, " +\
            f"indels: {numINDEL}, lowQual: {numLowQual}, multiAlt: {numMultiAlt}, good: {numGoodSNP}", file=sys.stderr)


vcf = VCF(snakemake.input['vcf'], snakemake.output[0], snakemake.params['qual_samtools'], 
          snakemake.params['heterzygote_cutoff'])
vcf.toNIEHS()


#if __name__ == '__main__':
#    with open("../20200429_BCFTOOLS/strain.order.vcf.txt") as s:
#        strains = s.read().strip().split()
#    # inputs and outputs
#    chromosome = list(range(1, 20)) + [ "X"]
#    for i in chromosome:
#        invcf = f"VCFs/combined.chr{i}.raw.vcf"
#        output = f"SNPs2/chr{i}.txt"
#        vcf = VCF(invcf, output, 50,  20, strains)
#        vcf.toNIEHS()
