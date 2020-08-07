import sys
import os
import pickle
import tqdm

codons = {"TTT":"F", "TTC":"F", "TCT":"S", "TCC":"S", "TAT":"Y", "TAC":"Y", "TGT":"C", "TGC":"C", 
          "TTA":"L", "TCA":"S", "TAA":"X", "TGA":"X", "TTG":"L", "TCG":"S", "TAG":"X", "TGG":"W", 
          "CTT":"L", "CTC":"L", "CCT":"P", "CCC":"P", "CAT":"H", "CAC":"H", "CGT":"R", "CGC":"R", 
          "CTA":"L", "CTG":"L", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGG":"R", "ATT":"I", 
          "ATC":"I", "ACT":"T", "ACC":"T", "AAT":"N", "AAC":"N", "AGT":"S", "AGC":"S", "ATA":"I", "ACA":"T", 
          "AAA":"K", "AGA":"R", "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R", "GTT":"V", "GTC":"V", "GCT":"A", 
          "GCC":"A", "GAT":"D", "GAC":"D", "GGT":"G", "GGC":"G", "GTA":"V", "GTG":"V", "GCA":"A", "GCG":"A", 
          "GAA":"E", "GAG":"E", "GGA":"G", "GGG":"G"}

f = "mm10/mm10_knownGene.txt"

input_chrome = sys.argv[2]

gene_locs = {x.rstrip().split('\t')[0]:"%s:%s-%s" %(x.rstrip().split('\t')[1], x.rstrip().split('\t')[3], x.rstrip().split('\t')[4]) for x in open(f)}
exon_len = {x.rstrip().split('\t')[0]:len(x.rstrip().split('\t')[8][:-1].split(',')) for x in open(f)}

annovar_output_dir = os.getcwd() 
by_chr = {}
by_trans = {}
trans2gene = {}

temp = {}
for line in open(os.path.join(annovar_output_dir, "%s.exonic_annotated" %(input_chrome))):
    line = line.rstrip().split('\t')
    loc = int(line[4]) 
    if loc not in temp:
        temp[loc] = {}
    for x in line[2].strip(',').split(','):
        if ':' not in x:
            continue
        x = x.split(':')
        t_name = x[1]
        if t_name not in trans2gene:
            trans2gene[t_name] = x[0]
        if t_name not in by_trans:
            by_trans[t_name] = {}
        c_loc = int(x[3][3:-1])
        by_trans[t_name][c_loc] = x
        by_trans[t_name][c_loc].append(input_chrome.split('r')[-1])
        by_trans[t_name][c_loc].append(loc)
        if t_name not in temp[loc]:
            temp[loc][t_name] = []
        temp[loc][t_name].append(x[2:]) #exon number, coding change, aa change, before codon, after codon
by_chr[input_chrome] = temp


need_to_check_regions_by_chr = {}
no_check_by_chr ={}

for i in by_trans:
    if len(by_trans[i]) == 1: #only 1 
        first = list(by_trans[i].keys())[0]
        chrome = by_trans[i][first][-2]
        if chrome not in no_check_by_chr:
            no_check_by_chr[chrome] = {}
        if by_trans[i][first][-1] not in no_check_by_chr[chrome]:
            no_check_by_chr[chrome][by_trans[i][first][-1]] = []
        if by_trans[i][first] not in no_check_by_chr[chrome][by_trans[i][first][-1]]:
            no_check_by_chr[chrome][by_trans[i][first][-1]].append(by_trans[i][first])
        continue
    keys = sorted(list(by_trans[i].keys()))
    check_for_this_trans = set()
    for j in range(len(keys) - 2):
        first = keys[j]
        second = keys[j+1]
        third = keys[j+2]
        if by_trans[i][first][-1] in check_for_this_trans:
            continue
        if first + 1 == second: ## consecutive
            first_name = "%s_%s_%s" %(by_trans[i][first][2], by_trans[i][first][4][:-1], by_trans[i][first][5])
            second_name = "%s_%s_%s" %(by_trans[i][second][2], by_trans[i][second][4][:-1], by_trans[i][second][5])
            third_name = "%s_%s_%s" %(by_trans[i][third][2], by_trans[i][third][4][:-1], by_trans[i][third][5])
            if first_name == second_name: #in the same codon
                if second + 1 == third and second_name == third_name: # triple
                    new_aa = {by_trans[i][first][-1]:[by_trans[i][first][3][2], by_trans[i][first][3][-1]], 
                              by_trans[i][second][-1]:[by_trans[i][second][3][2], by_trans[i][second][3][-1]], 
                              by_trans[i][third][-1]:[by_trans[i][third][3][2], by_trans[i][third][3][-1]]}
                    order = [by_trans[i][first][-1], by_trans[i][second][-1], by_trans[i][third][-1]]
                    check_for_this_trans.add(by_trans[i][first][-1])
                    check_for_this_trans.add(by_trans[i][second][-1])
                    check_for_this_trans.add(by_trans[i][third][-1])
                    if by_trans[i][first][-1] > by_trans[i][second][-1]:
                        order.append(-1)
                    else:
                        order.append(1)
                else:
                    for l in range(len(by_trans[i][first][5]) - 1):
                        if by_trans[i][first][5][l] == by_trans[i][first][3][2] and by_trans[i][first][6][l] == by_trans[i][first][3][-1]:
                            if by_trans[i][second][5][l+1] == by_trans[i][second][3][2] and by_trans[i][second][6][l+1] == by_trans[i][second][3][-1]:
                                # correct phase
                                break
                                
                    if l == 0:
                        new_aa = {by_trans[i][first][-1]:[by_trans[i][first][3][2], by_trans[i][first][3][-1]], 
                                  by_trans[i][second][-1]:[by_trans[i][second][3][2], by_trans[i][second][3][-1]], 
                                  "unchanged": [by_trans[i][first][5][-1], by_trans[i][first][5][-1]]}
                        order = [by_trans[i][first][-1], by_trans[i][second][-1], "unchanged"]
                        check_for_this_trans.add(by_trans[i][first][-1])
                        check_for_this_trans.add(by_trans[i][second][-1])
                        if by_trans[i][first][-1] > by_trans[i][second][-1]:
                            check_for_this_trans.add(by_trans[i][second][-1] - 1)
                            order[-1] = by_trans[i][second][-1] - 1
                            order.append(-1)
                        else:
                            check_for_this_trans.add(by_trans[i][second][-1] + 1)
                            order[-1] = by_trans[i][second][-1] + 1
                            order.append(1)
                    elif l == 1:
                        new_aa = {"unchanged": [by_trans[i][first][5][0], by_trans[i][first][5][0]],
                                  by_trans[i][first][-1]:[by_trans[i][first][3][2], by_trans[i][first][3][-1]], 
                                  by_trans[i][second][-1]:[by_trans[i][second][3][2], by_trans[i][second][3][-1]]} 
                        order = ["unchanged", by_trans[i][first][-1], by_trans[i][second][-1]] 
                        check_for_this_trans.add(by_trans[i][first][-1])
                        check_for_this_trans.add(by_trans[i][second][-1])
                        if by_trans[i][first][-1] > by_trans[i][second][-1]:
                            check_for_this_trans.add(by_trans[i][first][-1] + 1)
                            order[0] = by_trans[i][first][-1] + 1
                            order.append(-1)
                        else:
                            check_for_this_trans.add(by_trans[i][first][-1] - 1)
                            order[0] = by_trans[i][first][-1] - 1
                            order.append(1)
                    else: #error
                        print("ERROR")
                chrome = by_trans[i][first][-2]
                e_loc = by_trans[i][first][2]
                if chrome not in need_to_check_regions_by_chr:
                    need_to_check_regions_by_chr[chrome] = {}
                for l in new_aa:
                    if l == "unchanged":
                        continue
                    if l not in need_to_check_regions_by_chr[chrome]:
                        need_to_check_regions_by_chr[chrome][l] = {}
                    need_to_check_regions_by_chr[chrome][l][i] = [new_aa, order, e_loc]
            else: # just keep 
                chrome = by_trans[i][first][-2]
                if chrome not in no_check_by_chr:
                    no_check_by_chr[chrome] = {}
                if by_trans[i][first][-1] not in no_check_by_chr[chrome]:
                    no_check_by_chr[chrome][by_trans[i][first][-1]] = []
                if by_trans[i][first] not in no_check_by_chr[chrome][by_trans[i][first][-1]]:
                    no_check_by_chr[chrome][by_trans[i][first][-1]].append(by_trans[i][first])
        else:
            chrome = by_trans[i][first][-2]
            if chrome not in no_check_by_chr:
                no_check_by_chr[chrome] = {}
            if by_trans[i][first][-1] not in no_check_by_chr[chrome]:
                no_check_by_chr[chrome][by_trans[i][first][-1]] = []
            if by_trans[i][first] not in no_check_by_chr[chrome][by_trans[i][first][-1]]:
                no_check_by_chr[chrome][by_trans[i][first][-1]].append(by_trans[i][first])
                
    first = keys[-1]
    chrome = by_trans[i][first][-2]
    if by_trans[i][first][-1] not in check_for_this_trans:
        if chrome not in no_check_by_chr:
            no_check_by_chr[chrome] = {}
        if by_trans[i][first][-1] not in no_check_by_chr[chrome]:
            no_check_by_chr[chrome][by_trans[i][first][-1]] = []
        if by_trans[i][first] not in no_check_by_chr[chrome][by_trans[i][first][-1]]:
            no_check_by_chr[chrome][by_trans[i][first][-1]].append(by_trans[i][first])
        
    first = keys[-2]
    chrome = by_trans[i][first][-2]
    if by_trans[i][first][-1] not in check_for_this_trans:
        if chrome not in no_check_by_chr:
            no_check_by_chr[chrome] = {}
        if by_trans[i][first][-1] not in no_check_by_chr[chrome]:
            no_check_by_chr[chrome][by_trans[i][first][-1]] = []
        if by_trans[i][first] not in no_check_by_chr[chrome][by_trans[i][first][-1]]:
            no_check_by_chr[chrome][by_trans[i][first][-1]].append(by_trans[i][first])


## LOAD all SNPS
snp_file = os.path.join(sys.argv[1], "%s.txt" %input_chrome) # "/cluster/u/byoo1/jobs/HBCGM/peltz_code/HBCGM/PELTZ_20190301/SNPS"
real_order = []
by_chr_snp = {}
for line in tqdm.tqdm(open(snp_file)):
    if not real_order:
        real_order = line.rstrip().split('\t')
        continue
    if line.startswith("C57BL/6J"):
        continue
    line = line.rstrip().split('\t')
    if line[1] not in by_chr_snp:
        by_chr_snp[line[1]] = {}
    loc = int(line[2])
    if loc not in by_chr_snp[line[1]]:
        by_chr_snp[line[1]][loc] = line[-1]

other_variants = {}
for line in open("gene_annotation_%s.txt" %input_chrome):
    if line.startswith("{"):
        continue
    line = line.rstrip().split('\t')
    if len(line) < 3:
        continue
    chrome = line[0].split('_')[1]
    loc = int(line[0].split('_')[2])
    if chrome not in other_variants:
        other_variants[chrome] = {}
    other_variants[chrome][loc] = {}
    for case in range(1, len(line), 2):
        g_name = line[case]
        sym = line[case + 1]
        m = '!'.join([s for s in sym.split('!') if "<->" not in s and "CODING" not in s])
        if m:
            other_variants[chrome][loc][g_name] = m

def reverse( c ):
    c = c.upper()
    r = ""
    for i in c:
        if i == "A":
            r += "T"
        elif i == "T":
            r +=  "A"
        elif i == "C":
            r += "G"
        else:
            r += "C"
    return r

trans_dir = {x.split('\t')[0]:x.split('\t')[2] for x in open("mm10/mm10_knownGene.txt")}
gene2trans = {x.split('\t')[4]:x.split('\t')[0] for x in open("mm10/mm10_kgXref.txt")}
trans2gene = {x.split('\t')[0]:x.split('\t')[4] for x in open("mm10/mm10_kgXref.txt")}

AA_by_strains = {}

#codon tag if x %3 == 0 == last, == 1 == first, == 2 middle e.g. c.G771A , 771 % 3 0 so last GAG:GAA
for chrome in tqdm.tqdm(by_chr_snp):
    AA_by_strains[chrome] = {}
    temp = need_to_check_regions_by_chr[chrome] #multiple in codon
    temp2 = no_check_by_chr[chrome] #in exonic
    temp3 = other_variants[chrome] # other variants 
    checked = set()
    for snp in by_chr_snp[chrome]:

        length_strains = len(by_chr_snp[chrome][snp])
        AA_by_strains[chrome][snp] = {}
        if snp in temp: 

            for trans in temp[snp]:
                use = []
                if 'unchanged' in temp[snp][trans][0]:
                    use = ''.join([temp[snp][trans][0]['unchanged'][0] for x in range(length_strains)])
                
                s1 = by_chr_snp[chrome][temp[snp][trans][1][0]] if temp[snp][trans][1][0] in by_chr_snp[chrome] else use
                s2 = by_chr_snp[chrome][temp[snp][trans][1][1]] if temp[snp][trans][1][1] in by_chr_snp[chrome] else use
                s3 = by_chr_snp[chrome][temp[snp][trans][1][2]] if temp[snp][trans][1][2] in by_chr_snp[chrome] else use
                l1 = temp[snp][trans][1][0]
                l2 = temp[snp][trans][1][1]
                l3 = temp[snp][trans][1][2]
                if temp[snp][trans][1][-1] == -1: 
                    if l1 in by_chr_snp[chrome]:
                        s1 = reverse(s1)
                    if l2 in by_chr_snp[chrome]:
                        s2 = reverse(s2)
                    if l3 in by_chr_snp[chrome]:
                        s3 = reverse(s3)
                c = []
                
                for i in range(length_strains):
                    if s1[i] == "?" or s2[i] == "?" or s3[i] == "?":
                        c.append("?")
                    else:
                        c.append("%s%s%s" %(s1[i], s2[i], s3[i]))
                if "%s*%s" %(trans, l1) not in checked:
                    if l1 not in AA_by_strains[chrome]:
                        AA_by_strains[chrome][l1] = {}
                    if trans not in AA_by_strains[chrome][l1]:
                        AA_by_strains[chrome][l1][trans] = []
                    AA_by_strains[chrome][l1][trans].append(c)
                    
                if "%s*%s" %(trans, l2) not in checked:
                    if l2 not in AA_by_strains[chrome]:
                        AA_by_strains[chrome][l2] = {}
                    if trans not in AA_by_strains[chrome][l2]:
                        AA_by_strains[chrome][l2][trans] = []
                    AA_by_strains[chrome][l2][trans].append(c)
                    
                if "%s*%s" %(trans, l3) not in checked: 
                    if l3 not in AA_by_strains[chrome]:
                        AA_by_strains[chrome][l3] = {}
                    if trans not in AA_by_strains[chrome][l3]:
                        AA_by_strains[chrome][l3][trans] = []
                    AA_by_strains[chrome][l3][trans].append(c)
                checked.add("%s*%s" %(trans, l1))
                checked.add("%s*%s" %(trans, l2))
                checked.add("%s*%s" %(trans, l3))

        if snp in temp2:
            for trans in temp2[snp]:
                d = trans_dir[trans[1]]
                r = trans[3][2] if d == "+" else reverse(trans[3][2])
                a = trans[3][-1] if d == "+" else reverse(trans[3][-1])
                ref = trans[5]
                alt = trans[6]
                c = []
                l = temp2[snp][0][-1]
                s = by_chr_snp[chrome][l]
                for i in range(length_strains):
                    if s[i] == "?":
                        c.append("?")
                    else:
                        if s[i] == r:
                            c.append(ref)
                        else:
                            c.append(alt)
                if l not in AA_by_strains[chrome]:
                    AA_by_strains[chrome][l] = {}
                if trans[1] in AA_by_strains[chrome][l]:
                    print(c, l, AA_by_strains[chrome][l][trans[1]])
                if trans[1] not in AA_by_strains[chrome][l]:
                    AA_by_strains[chrome][l][trans[1]] = []
                AA_by_strains[chrome][l][trans[1]].append(c)
        if snp in temp3:
            if snp not in AA_by_strains[chrome]:
                AA_by_strains[chrome][snp] = {}
            for trans in temp3[snp]:
                if trans not in AA_by_strains[chrome][snp]:
                    AA_by_strains[chrome][snp][trans] = []
                AA_by_strains[chrome][snp][trans].append(temp3[snp][trans])

pickle.dump(AA_by_strains, open("AA_by_strains_%s.pkl" %input_chrome, 'wb'))


