//
// Created by Zhuoqing Fang on 2/1/22.
//

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include "dynum.h"
#include "ColumnReader.h"

struct Key {
    std::string chrom;
    int start;
    int end;
    Key(std::string key);
};

struct VEPSummary 
{
    std::string idx; //0
    std::string location; //1 /// for snp, its 'chr:pos', for indel,sv, its 'chr:start-end'
    std::string allel;// 2, for snp, indel -> dna seqence; sv -> "deletion", "duplicatoin" ...
    std::string geneid; //3 ensembl_gene_id
    std::string transxid; //4 ensembl_transcript_id
    std::string transtype; //5 -
    Dynum<std::string> consequence; //6 -
    std::string protein_position; //9  -
    Dynum<std::string> amino_acids; //10 -
    Dynum<std::string> codons; //11 -
    Dynum<std::string> dbsnpid; //12 -
    Dynum<std::string> samples; // sample name 13
    //std::string samples;
    std::string zygosity; // zygosity 14 HOM, HET
    Dynum<std::string> impact; //15
    std::string variant_class; //19 SO variant class: "SNV", "indel", deletion, duplication, inversion, ect
    std::string symbol; //20
    std::string biotype; //23, protein_coding, lincRNA, etc
    Dynum<std::string> HGVSc; //31 HGVSp strings
    Dynum<std::string> HGVSp; //32 HGVSp strings
    std::string chrom;
    int start;
    int end;

   VEPSummary(std::string uploaded_variant, std::string loc, std::string seq, std::string gene,
                std::string transcript, std::string feature_type, std::string csq, std::string aa_pos, std::string aa,
                std::string codon, std::string existing_variation, std::string ind, std::string zyg, std::string imp,
                std::string var_class, std::string gene_name, std::string bt, std::string hgvsc, std::string hgvsp);
};    




class VarirantEeffectPredictor
{
private:
    int numtoks;
    bool hasIND; // whether contain samples (IND) column in vep output
    std::unordered_map<std::string, int> PRIOR;
    std::unordered_map<std::string, std::string> CODONs;
    std::unordered_map<std::string, std::string> CSQs;
    Dynum<std::string> strainAbbrevs;
    std::string _line;
    // variant loc -> transxid -> record, orded by variant loc
    std::unordered_map<std::string, std::unordered_map<std::string, VEPSummary* >> geneCodingMap; 
    std::vector<VEPSummary*> data;
    std::unordered_map<std::string, unsigned> columns; // header colum name to index
    std::vector<std::string> keys; // output order sorted by key


    void lowercase(std::string & str);
    void upcase(std::string & str);
    void readHeader(char *inFileName, char *delemiter);
    bool compareKey(Key key1, Key key2);
    std::string codonChange(VEPSummary * pRecord); // reformat codon chagne strings e.g. CTA/L<->CTG/L
    std::string set_key(std::string location, std::string var_class);

public:
    VarirantEeffectPredictor(char* inVEPName, char* inStrainName);
    ~VarirantEeffectPredictor();
    void readVEP(char* inVEPName, char* delemiter, char* varType);
    void writeVEPImpact(char* outFileName);
    void writeVEPCsq(char* outFileName);
    bool getline(std::string &);
};