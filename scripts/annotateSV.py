import sys
import pandas as pd

from typing import Iterable, Tuple, List, Union, Optional, Dict

codons = {"TTT":"F", "TTC":"F", "TCT":"S", "TCC":"S", "TAT":"Y", "TAC":"Y", "TGT":"C", "TGC":"C", 
          "TTA":"L", "TCA":"S", "TAA":"X", "TGA":"X", "TTG":"L", "TCG":"S", "TAG":"X", "TGG":"W", 
          "CTT":"L", "CTC":"L", "CCT":"P", "CCC":"P", "CAT":"H", "CAC":"H", "CGT":"R", "CGC":"R", 
          "CTA":"L", "CTG":"L", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGG":"R", "ATT":"I", 
          "ATC":"I", "ACT":"T", "ACC":"T", "AAT":"N", "AAC":"N", "AGT":"S", "AGC":"S", "ATA":"I", "ACA":"T", 
          "AAA":"K", "AGA":"R", "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R", "GTT":"V", "GTC":"V", "GCT":"A", 
          "GCC":"A", "GAT":"D", "GAC":"D", "GGT":"G", "GGC":"G", "GTA":"V", "GTG":"V", "GCA":"A", "GCG":"A", 
          "GAA":"E", "GAG":"E", "GGA":"G", "GGG":"G"}


CSQs = {
        "transcript_ablation":"HIGH",
        "splice_acceptor_variant":"HIGH",
        "splice_donor_variant":"HIGH",
        "stop_gained":"HIGH",
        "frameshift_variant":"HIGH",
        "stop_lost":"HIGH",
        "start_lost":"HIGH",
        "transcript_amplification":"HIGH",
        "inframe_insertion":"MODERATE",
        "inframe_deletion":"MODERATE",
        "missense_variant":"MODERATE",
        "protein_altering_variant":"MODERATE",
        "splice_region_variant":"LOW",
        "incomplete_terminal_codon_variant":"LOW",
        "start_retained_variant":"LOW",
        "stop_retained_variant":"LOW",
        "synonymous_variant":"LOW",
        "coding_sequence_variant":"MODIFIER",
        "mature_miRNA_variant":"MODIFIER",
        "5_prime_UTR_variant":"MODIFIER",
        "3_prime_UTR_variant":"MODIFIER",
        "non_coding_transcript_exon_variant":"MODIFIER",
        "intron_variant":"MODIFIER",
        "NMD_transcript_variant":"MODIFIER",
        "non_coding_transcript_variant":"MODIFIER",
        "upstream_gene_variant":"MODIFIER",
        "downstream_gene_variant":"MODIFIER",
        "TFBS_ablation":"MODIFIER",
        "TFBS_amplification":"MODIFIER",
        "TF_binding_site_variant":"MODIFIER",
        "regulatory_region_ablation":"MODERATE",
        "regulatory_region_amplification":"MODIFIER",
        "feature_elongation":"MODIFIER",
        "regulatory_region_variant":"MODIFIER",
        "feature_truncation":"MODIFIER",
        "intergenic_variant":"MODIFIER",
        }


def unique(seq: Iterable) -> List[str]:
    """Remove duplicates from a list in Python while preserving order.
    :param seq: a python list object.
    :return: a list without duplicates while preserving order.
    """
    seen = set()
    seen_add = seen.add

    return [x for x in seq if x not in seen and not seen_add(x)]

def info2dict(info):
    INFO = info.split(";")
    INFO2 = dict(item.split('=')  for item in INFO if len(item.split('='))==2)
    return INFO2
def vcf2dict(invcf: str):
    gt_dict = {}
    with open(invcf, 'r') as filtered:
        for line in filtered:
            if line.startswith("#"): 
                continue  
            record = line.strip().split("\t")
            ### filtering
            INFO = info2dict(record[7])
            CHROM = record[0]
            ID = record[2]
            svtype = INFO['SVTYPE']
            sv_size = int(float(INFO['SVLEN']))
            start = int(record[1]) # bed coordinate
            end = int(INFO['END'])
            gt_dict[ID] =f"SV_{CHROM}_{start}_{end}_{svtype}"

    return gt_dict




def get_annotation(vcf:str, vep:str, out: str, *args, **kwargs):
    header = ""
    with open(vep,'r') as v:
        for line in v:
            if line.startswith("#Uploaded_variation"):
                header = line
                break
    vep = pd.read_table(vep, comment="#", header=None, dtype=str)
    vep.columns = header.strip("#\n").split("\t")
    output =  open(out, 'w') 
    print("Read VCF")
    vcf_dict = vcf2dict(vcf)
    vep = vep.set_index('Uploaded_variation')
    print("Generate Annotation for SVs")
    for k, v in vcf_dict.items():
        if k not in vep.index:
            output.write(v + "\n")
            continue
        anno = vep.loc[k,['Feature_type', 'SYMBOL','Consequence','Protein_position']]
        if anno.empty:
            continue
        if isinstance(anno, pd.DataFrame):
             anno = anno.drop_duplicates()
        else:
            anno = anno.to_frame().T
        #anno['Consequence'] = anno['Consequence'].str.spit(",").str[0]
        anno = anno[anno['SYMBOL'] != "-"]
        if anno.empty:
            continue
        csq = anno['SYMBOL'] +"\t" + anno['Consequence'] + ":"+anno['Protein_position']
        csq_out = "\t".join(csq.to_list())
        output.write(v + "\t" + csq_out + "\n")


    output.close()


def to_haplomap(annotation: str, outfile: str):
    output = open(outfile,'w')
    with open(annotation, 'r') as anno:
        for line in anno:
            ann = line.strip("\n").split("\t")
            for i in range(2, len(ann), 2):
                a = ann[i].split(":")[0].split(",")[0]
                if a in CSQs:
                    bit = CSQs[a]
                    ann[i] = bit
            outline = "\t".join(ann)
            output.write(outline+"\n")
    output.close()

if __name__ == "__main__":
    # inVCF = "/data/bases/fangzq/20200815_SV/svtools/merged.sv.pruned.filtered.vcf"
    # inVEP = "/data/bases/fangzq/20200815_SV/svtools/merged.sv.pruned.filtered.pass.vep.txt"
    # outTXT = "/data/bases/fangzq/20200815_SV/svtools/SV_annotation.txt"
    # outBlockAnno = "/data/bases/fangzq/20200815_SV/svtools/SV_annotation_eblocks.txt" ## for eblock input "-g"
    inVCF = sys.argv[1]
    inVEP = sys.argv[2]
    outTXT = sys.argv[3]
    outBlockAnno = sys.argv[4]

    # get_annotation(inVCF, inVEP, outTXT)
    to_haplomap(outTXT, outBlockAnno)
    print("Done")


        
        


