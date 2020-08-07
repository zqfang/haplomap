#python check.py mm10_fastafile mm10_knownGene.txt new_fastafile
import os
import shutil
import sys

c_start_codon = "ATG"
c_stop_codon = ["TAA", "TGA", "TAG"]

fa = open(os.path.join(sys.argv[1], "mm10_knownGeneMrna.fa" ))
g = open(os.path.join(sys.argv[1], "mm10_knownGene.txt"))

fa_dict = {}
curr = ""
for line in fa:
    if line.startswith('>'):
        curr = line.split()[0].strip('>')
        continue   
    if curr not in fa_dict: 
        fa_dict[curr] = line.rstrip()
    else:
        fa_dict[curr] += line.rstrip()

save_out = open("saved_transcript.txt", 'w') 
toss_out = open("tossed_transcript.txt", 'w') 
saved = set()

for line in g:
    name,Chr,dbstrand,txstart,txend,cdsstart,cdsend,exoncount,exonstart,exonend = line.rstrip().split('\t')[:10]
    if name not in fa_dict:
        continue
    use = fa_dict[name]

    start_loc = int(txstart) + 1
    end_loc = int(txend)
    cds_start = int(cdsstart) + 1
    cds_end = int(cdsend)
    
    if cds_start == (cds_end + 1): #no coding sequence
        continue

    exons = []
    seq_length = 0
    for i in range(len(exonstart.split(',')) -1):
        exons.append([int(exonstart.split(',')[i]), int(exonend.split(',')[i])])
        seq_length += (int(exonend.split(',')[i]) - int(exonstart.split(',')[i]))
    if dbstrand == "+":
        curr = [(x[0] + 1)  for x in exons if x[0] <= cds_start and x[1] >= cds_start ][0]
        pre_length = sum([x[1]-x[0] for x in exons if x[1] < curr])
        start = cds_start - curr + pre_length
        start_codon = use[start:start+3]
        curr = [x[0] for x in exons if x[0] <= cds_end and x[1] >= cds_end ][0]
        pre_length = sum([x[1]-x[0] for x in exons if x[1] < curr])
        end = pre_length + (cds_end - curr)
        end_codon = use[end-3:end]
    else:
        # the fasta is already reversed complemented 
        curr = [x[1]  for x in exons if x[0] <= cds_end and x[1] >= cds_end ][0]
        pre_length = sum([x[1]-x[0] for x in exons if x[1] > curr])
        start = curr - cds_end + pre_length
        start_codon = use[start:start+3]
        curr = [(x[1] + 1) for x in exons if x[0] <= cds_start and x[1] >= cds_start ][0]
        pre_length = sum([x[1]-x[0] for x in exons if x[0] > curr])
        end = pre_length + ( curr - cds_start) # - curr)
        end_codon = use[end-3:end]
    
    if start_codon == c_start_codon and end_codon in c_stop_codon:
        save_out.write("%s\n" %(name))
        saved.add(name)
    else:
        toss_out.write("%s\n" %(name))
save_out.close()
toss_out.close()


exons = {}
g = open(os.path.join(sys.argv[1], "mm10_knownGene.txt"))
for line in g:
    name,Chr,dbstrand,txstart,txend,cdsstart,cdsend,exoncount,exonstart,exonend = line.rstrip().split('\t')[:10]
    exons[name] = "Comment: this sequence (leftmost exon at %s:%s) is generated" %(Chr, exonstart.split(',')[0])


curr = ""
l = ""
out = open("new_fasta", 'w')

for trans in saved:
    out.write(">%s\t%s\n" %(trans, exons[trans]))
    out.write("%s\n" %(fa_dict[trans]))
out.close()

out = open("new_knownGene", 'w')
g = open(os.path.join(sys.argv[1], "mm10_knownGene.txt"))
for line in g:
    name,Chr,dbstrand,txstart,txend,cdsstart,cdsend,exoncount,exonstart,exonend = line.rstrip().split('\t')[:10] 
    if name not in saved:
        continue
    out.write(line)
out.close()

if not os.path.isdir("mm10"):
    os.mkdir("mm10")
shutil.move("new_fasta", "mm10/mm10_knownGeneMrna.fa")
shutil.move("new_knownGene", "mm10/mm10_knownGene.txt")




