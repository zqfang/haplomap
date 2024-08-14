
#include "constants.h"

std::unordered_map<std::string, std::string> CODONs = {{"TTT", "F"}, {"TTC", "F"}, {"TCT", "S"}, {"TCC", "S"}, {"TAT", "Y"}, {"TAC", "Y"}, {"TGT", "C"}, {"TGC", "C"}, 
          {"TTA", "L"}, {"TCA", "S"}, {"TAA", "X"}, {"TGA", "X"}, {"TTG", "L"}, {"TCG", "S"}, {"TAG", "X"}, {"TGG", "W"}, 
          {"CTT", "L"}, {"CTC", "L"}, {"CCT", "P"}, {"CCC", "P"}, {"CAT", "H"}, {"CAC", "H"}, {"CGT", "R"}, {"CGC", "R"}, 
          {"CTA", "L"}, {"CTG", "L"}, {"CCA", "P"}, {"CCG", "P"}, {"CAA", "Q"}, {"CAG", "Q"}, {"CGA", "R"}, {"CGG", "R"}, {"ATT", "I"}, 
          {"ATC", "I"}, {"ACT", "T"}, {"ACC", "T"}, {"AAT", "N"}, {"AAC", "N"}, {"AGT", "S"}, {"AGC", "S"}, {"ATA", "I"}, {"ACA", "T"}, 
          {"AAA", "K"}, {"AGA", "R"}, {"ATG", "M"}, {"ACG", "T"}, {"AAG", "K"}, {"AGG", "R"}, {"GTT", "V"}, {"GTC", "V"}, {"GCT", "A"}, 
          {"GCC", "A"}, {"GAT", "D"}, {"GAC", "D"}, {"GGT", "G"}, {"GGC", "G"}, {"GTA", "V"}, {"GTG", "V"}, {"GCA", "A"}, {"GCG", "A"}, 
          {"GAA", "E"}, {"GAG", "E"}, {"GGA", "G"}, {"GGG", "G"}};

 
/// see the variant impact coding here: 
/// https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
/// last update: 2024-08-15
std::unordered_map<std::string, int> PRIOR = {{"HIGH", 2}, {"MODERATE", 1}, {"LOW", 0}, {"MODIFIER", -1}};
std::unordered_map<std::string, int> CSQs = {
    {"transcript_ablation", 2},
    {"splice_acceptor_variant", 2},
    {"splice_donor_variant", 2},
    {"stop_gained", 2},
    {"frameshift_variant", 2},
    {"stop_lost", 2},
    {"start_lost", 2},
    {"transcript_amplification", 2},
    {"feature_elongation", 2},
    {"feature_truncation", 2},
    {"inframe_insertion", 1},
    {"inframe_deletion", 1},
    {"missense_variant", 1},
    {"protein_altering_variant", 1},
    {"splice_donor_5th_base_variant", 0},
    {"splice_region_variant", 0},
    {"splice_donor_region_variant", 0},
    {"splice_polypyrimidine_tract_variant", 0},
    {"incomplete_terminal_codon_variant", 0},
    {"start_retained_variant", 0},
    {"stop_retained_variant", 0},
    {"synonymous_variant", 0},
    {"coding_sequence_variant", -1},
    {"mature_miRNA_variant", -1},
    {"5_prime_UTR_variant", -1},
    {"3_prime_UTR_variant", -1},
    {"non_coding_transcript_exon_variant", -1},
    {"intron_variant", -1},
    {"NMD_transcript_variant", -1},
    {"non_coding_transcript_variant", -1},
    {"coding_transcript_variant", -1},
    {"upstream_gene_variant", -1},
    {"downstream_gene_variant", -1},
    {"TFBS_ablation", -1},
    {"TFBS_amplification", -1},
    {"TF_binding_site_variant", -1},
    {"regulatory_region_ablation", -1},
    {"regulatory_region_amplification", -1},
    {"regulatory_region_variant", -1},
    {"intergenic_variant", -1},
    {"sequence_variant", -1},
};