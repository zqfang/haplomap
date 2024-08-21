//
// Created by Zhuoqing Fang on 2/22/21.
//
#ifndef VCF_H
#define VCF_H

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include "utils.h"
#include "dynum.h"


struct VCFOptions
{
    char *inputFileName;
    char *outputFileName;
    char *sampleFileName;
    char *variantType; 
    bool plink; 
    float phredLikelihoodDifference; // PL_other - PL_homo > 20
    float qual;  // QUAL field
    float mappingQuality; // MQ
    float strandBiasPhredPval; // SP (BCFtools) or FS (GATK)
    float readPositionBias; //RPB, or endDistanceBias (ReadPosRankSum)
    float baseQualityBias; //BQB
    float variantDistanceBias; //VDB
    /// allele depth of coverage (AD)
    /// variants < AD set to be LowQual
    float alleleDepth;
    /// see here https://samtools.github.io/bcftools/howtos/variant-calling.html
    /// ratio = %MAX(AD) / %MAX(DP). variants < ratio set to be LowQual
    float ratio; 
    int gapWindow;
    int snpGap;
    bool homozygous; 
    bool verbose;

    // constructor
    VCFOptions() : inputFileName(nullptr), outputFileName(nullptr), 
                     sampleFileName(nullptr), variantType((char*)"SNV"),
                     plink(false), phredLikelihoodDifference(20.0), 
                     qual(50.0), mappingQuality(20.0), strandBiasPhredPval(50.0), 
                     readPositionBias(0.0001), baseQualityBias(0),variantDistanceBias(0), 
                     alleleDepth(3.0), ratio(0.1), gapWindow(20), snpGap(3), 
                     homozygous(true), verbose(false) {};
};

class Variant {

public:
    int TYPE; // snv: 1, indel: 2, sv: 3, unknown: 0
    std::string CHROM;
    int POS;
    std::string ID;
    std::string REF;
    std::string ALT;
    float QUAL;
    std::string FILTER;
    std::string INFO;
    std::unordered_map<std::string, std::vector<std::string>> FORMATS;
     
    Variant() = default;
    Variant(std::string record);
    std::unordered_map<std::string, std::string> getINFO();
    float maxAllelicDepth();
    float maxStrandBiasPhredScalePval(); 
    float maxDepth();

private:
    void parseFORMAT(const std::vector<std::string> & rec);
    bool isSNP();
    bool isSV(); // is structral variant ?
    bool isINDEL();
    float max(const std::string & fmt);  
};



class VCF 
{
public:
    VCF() = default;
    VCF(std::istream* input, std::shared_ptr<VCFOptions> opts);
    VCF(std::istream* input,  std::shared_ptr<VCFOptions> opts, std::vector<std::string> samples);
    ~VCF();
    void parseRecords();
    void parseHeader();
    void checkSamples();
    bool parseStructralVariant(std::string& alleles, std::vector<std::string>& alts, std::vector<int>& hasAlt);
    bool parseSNP(std::string& alleles, std::vector<std::string>& alts, std::vector<int>& hasAlt);
    void writeSNP(std::string& alleles);
    void writeTPED(std::string& alleles);
    void writeTFAM();
    void writeStructralVariant(std::string& alleles);

private:
    std::string line;
    Variant variant;
    std::string outHeader;
    Dynum<std::string> strains;
    std::vector<std::string> samples;
    // output file
    std::ofstream output;
    std::ofstream tped;
    std::ofstream tfam;
    // read input file
    std::istream* input;
    VCFOptions* opts;
    char getAllele(std::string& alleles, int i) { return alleles[i]; }
    void setAllele(std::string& alleles, int i, char a) { alleles[i] = a; }
};

#endif