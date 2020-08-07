import sys
import tqdm


regdoms = {}
for i in open("mm10_regdoms"):
    Chr,Start,End,Gene,V=i.rstrip().split('\t')
    Chr = Chr.split('r')[-1]
    Start = int(Start)
    End = int(End)
    if Chr not in regdoms:  
        regdoms[Chr] = {}
    if Start not in regdoms[Chr]:
        regdoms[Chr][Start] = {}
    if End not in regdoms[Chr][Start]:
        regdoms[Chr][Start][End] = []
    regdoms[Chr][Start][End].append(Gene)

def find_intergenic( Chr, Pos, regdoms=regdoms):
    return_list = []
    if Chr not in regdoms:
        return return_list
    keys = regdoms[Chr].keys()
    for Start in sorted(keys):
        if Start > Pos:
            continue
        for End in regdoms[Chr][Start]:
            if End >= Pos:
                return_list = return_list + regdoms[Chr][Start][End]
    return return_list

def get_change( gene, aachange ):
    changes = set()
    for i in aachange:
        i = i.split(':')
        if i[0] == gene:
            allele = ""
            protein = ""
            for j in i:
                if j.startswith("c."):
                    allele = [ j.split('.')[-1][0], j.split('.')[-1][-1] ]
                if j.startswith("p."):
                    protein = [ j.split('.')[-1][0], j.split('.')[-1][-1] ]
            print(aachange, gene, protein)
            if protein[0] == protein[1]:
                changes.add( "SYNONYMOUS_CODING")
            else:
                changes.add("%s/%s<->%s/%s" %(allele[0], protein[0], allele[1], protein[1] ))
    return changes

sig = set(["EXONIC", "SPLICE_SITE", "UTR3", "UTR5"])
use_file = open(sys.argv[1]).readlines()
for i in tqdm.tqdm(use_file):
    if i.startswith("Chr"):
        continue
    Chr,Start,End,Func,Gene,AAChange = i.rstrip().split('\t')
    name = "SNP_%s_%s" %(Chr, Start)
    intergenics = find_intergenic( Chr, int(Start) )
    Gene = Gene.split(';') 
    Func = Func.upper().replace("UTR5", "5PRIME_UTR").replace("UTR3", "3PRIME_UTR").replace("SPLICING", "SPLICE_SITE").split(';')
    AAChange = AAChange.split(',')
    to_print = {}
    for num in range(len(Gene)):
        function = Func[num]
        gene = Gene[num]
        if function == "EXONIC":
            change = get_change(gene, AAChange)
            if change:
                if gene not in to_print:
                    to_print[gene] = set()
                to_print[gene].update(change)
            else:
                if gene not in to_print:
                    to_print[gene] = set()
                to_print[gene].add("misc_genic")
        elif "UTR" in function or "INTRONIC" == function or "SPLICE_SITE" in function:
            if gene not in to_print:
                to_print[gene] = set()
            to_print[gene].add(function)
    if not (sig.intersection(Func)): #add intergenic
        for gene in intergenics:
            if gene not in Gene:
                if gene not in to_print:
                    to_print[gene] = set()
                to_print[gene].add( "intergenic")
    for c in to_print:
        name = "%s\t%s\t%s" %(name, c, "!".join(to_print[c]))
        print(name)


