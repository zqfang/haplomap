import os
import sys
import pickle

codons = {"TTT":"F", "TTC":"F", "TCT":"S", "TCC":"S", "TAT":"Y", "TAC":"Y", "TGT":"C", "TGC":"C", 
          "TTA":"L", "TCA":"S", "TAA":"X", "TGA":"X", "TTG":"L", "TCG":"S", "TAG":"X", "TGG":"W", 
          "CTT":"L", "CTC":"L", "CCT":"P", "CCC":"P", "CAT":"H", "CAC":"H", "CGT":"R", "CGC":"R", 
          "CTA":"L", "CTG":"L", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGG":"R", "ATT":"I", 
          "ATC":"I", "ACT":"T", "ACC":"T", "AAT":"N", "AAC":"N", "AGT":"S", "AGC":"S", "ATA":"I", "ACA":"T", 
          "AAA":"K", "AGA":"R", "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R", "GTT":"V", "GTC":"V", "GCT":"A", 
          "GCC":"A", "GAT":"D", "GAC":"D", "GGT":"G", "GGC":"G", "GTA":"V", "GTG":"V", "GCA":"A", "GCG":"A", 
          "GAA":"E", "GAG":"E", "GGA":"G", "GGG":"G"}


def get_annotation(strains, snpdb, annodb, kgxref, knowngene, ofile1, ofile2):


    # trans_dir = {x.split('\t')[0]:x.split('\t')[2] for x in open("mm10_knownGene.txt")}
    # gene2trans = {x.split('\t')[4]:x.split('\t')[0] for x in open("mm10_kgXref.txt")}
    # trans2gene = {x.split('\t')[0]:x.split('\t')[4] for x in open("mm10_kgXref.txt")}
    trans_dir = {x.split('\t')[0]:x.split('\t')[2] for x in open(knowngene)}
    gene2trans = {x.split('\t')[4]:x.split('\t')[0] for x in open(kgxref)}
    trans2gene = {x.split('\t')[0]:x.split('\t')[4] for x in open(kgxref)}    

    #AA_by_strains = pickle.load(open("AA_by_strains.pkl", 'rb'))
    AA_by_strains = pickle.load(open(annodb, 'rb'))
    # real_order = open("chr1.txt").readlines()[0].rstrip().split('\t') #replace with the snp file
    real_order = open(snpdb).readlines()[0].rstrip().split('\t')
    f = os.path.join(strains) # strain file name
    new_order = [x.rstrip().split('\t')[0] for x in open(f)]
    ref_index = real_order.index('C57BL/6J')
    by_case = {}

    for chrome in AA_by_strains:
        by_case[chrome] = {}
        for snp in AA_by_strains[chrome]:
            by_case[chrome][snp] = {}
            for trans in AA_by_strains[chrome][snp]:
                sym = AA_by_strains[chrome][snp][trans][0]
                if len(sym) == len(real_order): #AA change
                    ref = sym[ref_index]
                    change = set(sym[new_order.index(x)] for x in new_order if sym[new_order.index(x)] != "?")
                    updated_sym = ""
                    for c in change:
                        updated_sym += "%s/%s<->%s/%s!" %(ref, codons[ref], c, codons[c])
                    by_case[chrome][snp][trans] = updated_sym[:-1]
                else:
                    by_case[chrome][snp][trans] = sym

    out = open(ofile1, 'w') # write ensemble id anno
    out_name = open(ofile2, 'w') # write hgnc anno
    for chrome in by_case:
        for snp in by_case[chrome]:
            name = "SNP_%s_%s" %(chrome, snp)
            name2 = name
            for trans in by_case[chrome][snp]:
                name = "%s\t%s\t%s" %(name, trans, by_case[chrome][snp][trans]) 
                if trans in trans2gene:
                    name2 = "%s\t%s\t%s" %(name2, trans2gene[trans], by_case[chrome][snp][trans])
            out.write("%s\n" %name)
            out_name.write("%s\n" %name2)
    out.close()
    out_name.close()

# snpdb, annodb, kgxref, knowngene, ensemble, hgnc
strain2trait(snakemake.input['strains'], snakemake.input['snpdb'], snakemake.input['annodb'], 
             snakemake.input['kgxref'], snakemake.input['knowngene'], 
             snakemake.output['ensemble'], snakemake.output['hgnc'])