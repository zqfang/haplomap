import os
import sys
import numpy as np

codons = {"TTT":"F", "TTC":"F", "TCT":"S", "TCC":"S", "TAT":"Y", "TAC":"Y", "TGT":"C", "TGC":"C", 
          "TTA":"L", "TCA":"S", "TAA":"X", "TGA":"X", "TTG":"L", "TCG":"S", "TAG":"X", "TGG":"W", 
          "CTT":"L", "CTC":"L", "CCT":"P", "CCC":"P", "CAT":"H", "CAC":"H", "CGT":"R", "CGC":"R", 
          "CTA":"L", "CTG":"L", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGG":"R", "ATT":"I", 
          "ATC":"I", "ACT":"T", "ACC":"T", "AAT":"N", "AAC":"N", "AGT":"S", "AGC":"S", "ATA":"I", "ACA":"T", 
          "AAA":"K", "AGA":"R", "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R", "GTT":"V", "GTC":"V", "GCT":"A", 
          "GCC":"A", "GAT":"D", "GAC":"D", "GGT":"G", "GGC":"G", "GTA":"V", "GTG":"V", "GCA":"A", "GCG":"A", 
          "GAA":"E", "GAG":"E", "GGA":"G", "GGG":"G"}

temp = {}
for line in open(sys.argv[1]):
    line = line.rstrip().split('\t')
    loc = int(line[4]) 
    if loc not in temp:
        temp[loc] = {}
    for x in line[2].strip(',').split(','):
        if ':' not in x:
            continue
        x = x.split(':')
        t_name = x[1]
        if t_name not in temp[loc]:
            temp[loc][t_name] = []
        temp[loc][t_name].append(x[2:]) #exon number, coding change, aa change, before codon, after codon

def add_to_dict( d, loc, trans, info, new_aa ):
    if loc not in d:
        d[loc] = []
    d[loc].append("%s:%s:%s:%s%s" %(trans, info[0], info[1], info[2][:-1], codons[new_aa]))

l = list(temp.keys())
codon_merged = {}

#exonic function
#line1  nonsynonymous SNV   Lgsn:ENSMUST00000062560.13:exon5:c.A1174T:p.T392S,  1   31204012    31204012    A   T   het .
exonic_function = {}
for i in open(sys.argv[1]):
    i = i.rstrip().split('\t')
    if "frameshift" in i[1]:
        continue
    var = "%s^$^%s" %(i[3], i[4])
    exonic = i[2].strip(',')
    if var not in exonic_function:
        exonic_function[var] = []

    if int(i[4]) in codon_merged: 
        for j in codon_merged[int(i[4])]:
            name1 = ':'.join(j.split(':')[:-1])
            name2 = ':'.join(exonic.split(':')[1:4])
            if name1 == name2:
                exonic = "%s:%s" %(exonic.split(':')[0], j)
    exonic_function[var].append(exonic)

#variant function
variant_function = {}
for i in open(sys.argv[2]):
   #exonic  Lgsn    1   31204012    31204012    A   T   het .
    i = i.rstrip().split('\t')
    var = "%s^$^%s" %(i[2], i[3])
    if "(" in i[1]:
        gene = [i[1].split("(")[0].strip()]
    else:
        gene = [x.strip() for x in i[1].split(',')]
    func = i[0]
    if var not in variant_function:
        variant_function[var] = set()
    for g in gene:
        variant_function[var].add("%s^$^%s" %(g, func))

f = "/cluster/u/byoo1/jobs/HBCGM/peltz_code/HBCGM/PELTZ_20190301/SNPS/%s" %sys.argv[3]
for i in open(f):
    if i.startswith("C57BL/6J"):
        continue
    i = i.rstrip().split('\t')
    var = "%s^$^%s" %(i[1], i[2])
    if var not in variant_function:
        variant_function[var] = set()

#output 
out = open(sys.argv[4], 'w')
out.write("Chr\tStart\tEnd\tFunc.knownGene\tGene.knownGene\tAAChange.knownGene\n")

for i in variant_function:
    to_print = "%s\t%s\t%s\t" %(i.split("^$^")[0], i.split("^$^")[1], i.split("^$^")[1] )
    func = ""
    gene = ""
    if not variant_function[i]:
        func = ".,"
        gene = ".,"
    else:
        for j in variant_function[i]:
            func += "%s;" %j.split("^$^")[1]
            gene += "%s;" %j.split("^$^")[0]
    to_print += "%s\t%s\t" %(func[:-1], gene[:-1])
    if i in exonic_function:
        to_print += ','.join(exonic_function[i])
    else:
        to_print += '.'
    out.write("%s\n" %(to_print))


