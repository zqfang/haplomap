//
// Created by Zhuoqing Fang on 2/22/21.
//

#include "vcf.h"
#include "constants.h"


Variant::Variant(std::string record)
{
    std::vector<std::string> rec = split(trim(record), '\t');

    CHROM = rec[0];
    size_t tok0 = CHROM.find_first_not_of("chr");
    CHROM = CHROM.substr(tok0, CHROM.size() - tok0);
    POS = std::stoi(rec[1]);
    ID = rec[2];
    REF = rec[3];
    ALT = rec[4];
    QUAL = std::stof(rec[5]);
    FILTER = rec[6];
    INFO = rec[7];
    // format dict
    TYPE = 0;
    if (this->isSNP())
    {
        TYPE = 1; // snv
    } else if (this->isSV()) {
        TYPE = 3; // sv
    } else if ( this->isINDEL() ||  REF.size() != ALT.size() ) {
        TYPE = 2; // indel
    }
    parseFORMAT(rec);
    
}
std::unordered_map<std::string, std::string> Variant::getINFO() 
{
    std::unordered_map<std::string, std::string> infos;
    std::vector<std::string> info = split(this->INFO,';');
    std::vector<std::string> _info;
    for (auto &item: info) 
    {
        _info = split(item,'=');
        // add values to hashmap
        infos[_info[0]] = _info.size() == 2 ? _info[1]: "true";
    }    
    return infos;
}
float Variant::maxAllelicDepth()
{
    if (this->FORMATS.find("AD") == this->FORMATS.end())
        return 0; // if AD not found, then we should skip 
    
    float AD(0);
    for (auto & a: this->FORMATS["AD"])
    {
        // AD is a vector for gt (ref0, alt1, alt2)
        std::vector<std::string> ad = split(a, ',');
        for (auto it = ad.begin(); it != ad.end(); it++) 
        {
            if (*it != ".")
                AD = std::max(std::stof(*it), AD);
        }
    }           
    return AD;   
}

float Variant::maxDepth()
{
    return this->max("DP");
}
float Variant::maxStrandBiasPhredScalePval() 
{
    return this->max("SP");
}

void Variant::parseFORMAT(const std::vector<std::string> & rec) 
{
    std::vector<std::string> _fmt = split(rec[8], ':');
    std::vector<std::vector<std::string>> temp(_fmt.size(), std::vector<std::string>(rec.size()-9, "?"));
    
    for (size_t i = 9; i < rec.size(); i ++) 
    {
        std::vector<std::string> sampleValues = split(rec[i],':');
        for (size_t j=0; j < _fmt.size(); j++)
            // add values to hashmap
            /// FIXME: Not all sampleValues have equal length, push_back("-999") for temp
            if (j < sampleValues.size()) 
            {
                FORMATS[_fmt[j]].push_back(sampleValues[j]);
            }
            else
            {
                FORMATS[_fmt[j]].push_back("-999");
            }  
    }
    // std::cout<<"debug"<<std::endl;
    // for (size_t j=0; j < _fmt.size(); j++)
    //      FORMATS.insert(std::make_pair(_fmt[j], temp[j]));
    return;
}

bool Variant::isSNP()
{
    if (this->REF.size() > 1) return false;
    if (this->REF == "N") return false;
    // e.g. <DEL>
    if (this->ALT.find("<") != std::string::npos || this->ALT.find(">") != std::string::npos)
        return false;

    std::vector<std::string> alts = split(this->ALT, ',');
    for (auto &alt: alts)
    {
        if (alt.size() > 1)
            return false;
    }
    return true;
}

bool Variant::isSV()
{
    if (this->INFO.find("SVTYPE=") != std::string::npos)
        return true;
    return false;
}

bool Variant::isINDEL()
{
    if (this->INFO.find("INDEL") != std::string::npos)
        return true;
    return false;
}

float Variant::max(const std::string & fmt) 
{
    if (this->FORMATS.find(fmt) == this->FORMATS.end() )
        return 0; // if DP not found, then we should skip 
    
    float res(0);
    /// FIXME: if fmt is .
    for (auto & a: this->FORMATS[fmt])
    {   
        if (a != ".")
            res = std::max(std::stof(a), res);
    }
    return res;
}   


VCF::VCF(std::istream* input, std::shared_ptr<VCFOptions> opts):
outHeader("C57BL/6J")
{
    this->input = input;
    this->opts = opts.get();
    this->output.open(opts->outputFileName);      
    // output header
    parseHeader();
    output<<outHeader;
    for (int i=0; i < strains.size(); i++)
    {
        output<<"\t"<<strains.eltOf(i);
    }
    output<<std::endl;
    
    if (opts->plink)
    {
       this->tped.open(std::string(opts->outputFileName)+".tped");
       this->tfam.open(std::string(opts->outputFileName)+".tfam");
       this->writeTFAM();
    }

}
VCF::VCF(std::istream* input,  std::shared_ptr<VCFOptions> opts, std::vector<std::string> samples):
outHeader("C57BL/6J")
{
    this->input = input;
    this->opts = opts.get();
    this->samples = samples;
    this->output.open(opts->outputFileName);   
  
    parseHeader();
    // output header
    output<<outHeader;
    if (!samples.empty()) 
    {
        int refind = -1;
        for (size_t i = 0; i < samples.size(); i++)
        {
            if (samples[i] == outHeader) 
            {
                refind = i;
                continue;
            }
            output <<"\t"<<samples[i];
        }   
        // drop reference strain if in the samples 
        if (refind >= 0)
            samples.erase(samples.begin()+refind);
    } else
    {
        for (int i=0; i < strains.size(); i++)
            output<<"\t"<<strains.eltOf(i);
    }
    output<<std::endl;

    if (opts->plink)
    {
       std::string pout(opts->outputFileName);
       size_t pos = pout.find_last_of(".");
       this->tped.open(pout.substr(0, pos)+".tped");
       this->tfam.open(pout.substr(0, pos)+".tfam");
       this->writeTFAM();
    } 
}


VCF::~VCF() 
{
    // input.close();
    output.close();
    if (opts->plink)
    {
        tped.close();
        tfam.close();
    }
}

void VCF::parseHeader()
{
    // read vcf header
    while (std::getline(*input, line))
    {
        // starts with "##"
        if (line.find("##") == 0) 
        {
            continue;
        }
        if (line.find("#CHROM") == 0) 
        {
            std::vector<std::string> header = split(trim(line),'\t');
            //strains = std::vector<std::string>(header.begin()+9, header.end());
            for (auto itr = header.begin()+9; itr != header.end(); itr++)
                strains.addElementIfNew(*itr);
            break;
        }  
    }
    this->checkSamples();
}
void VCF::checkSamples() 
{
    if (!samples.empty()) 
    {
        if (strains.size() < (int) samples.size())
        {
            std::cerr<<"Input (--samples) size larger than vcf sample size !!!"
                        <<std::endl;
            std::exit(-1);
        }
        // check samples in strains
        for (auto &s: samples)
        {
            int inStrain = strains.hasIndex(s);
            if (inStrain == -1)
            {
                std::cerr<<"Input name: "<<s<<"  NOT in VCF samples"<<std::endl;
                std::exit(1);
            }
        }
    }
}


bool VCF::parseSNP(std::string & alleles, 
                    std::vector<std::string>& alts, std::vector<int>& hasAlt)
{
    if (variant.REF == "N")
        return false;

    std::unordered_map<std::string, std::string> INFO = variant.getINFO();
    if (INFO["MQ"].size() < 1 || std::stof(INFO["MQ"]) < opts->mappingQuality)
        return false;

    float strandBiasPhredPval;
   
    if (INFO.find("FS") != INFO.end()) 
    {
        strandBiasPhredPval = std::stof(INFO["FS"]);
    } 
    else if (variant.FORMATS.find("SP") != variant.FORMATS.end()) 
    {
        strandBiasPhredPval = variant.maxStrandBiasPhredScalePval();
    } 
    else 
    {
        std::cerr<<"Variant position "<<variant.CHROM<<":"<<variant.POS
                    <<"INFO Tag (SB or FS) is not Found ! Skip strand bias filtering"
                    <<std::endl;
        strandBiasPhredPval = 0;
    } 
    if (strandBiasPhredPval >  opts->strandBiasPhredPval)
        return false; // filter strand bias variant 

    float AD = variant.maxAllelicDepth();
    if (AD < opts->alleleDepth || (AD / variant.maxDepth()) < opts->ratio)
        return false;
    // start parsing samples
    for (int s=0; s < strains.size(); s++)
    {
        std::string gt = variant.FORMATS["GT"][s];
        if (gt == "./." || gt == ".|.")
            continue;  // allele == '?'
        // replace
        std::replace(gt.begin(), gt.end(), '|', '/');
        std::vector<std::string> gts = split(gt, '/');
        
        if (gts.size() != 2)
            continue;
        if (gts[0] != gts[1])
            continue; // allele == '?'

        // genotype 0/0, 0/1, 1/1, 1/2,2/2 
        std::vector<int> GTs;
        // string to int, vectorize
        std::transform(gts.begin(), gts.end(), std::back_inserter(GTs),
                        [](std::string &s) { return std::stoi(s); });
        unsigned int ind = (GTs[0] + 1) * (GTs[0] + 2) / 2 - 1;

        if (variant.FORMATS.find("PL") == variant.FORMATS.end())
        {
            std::cerr<<"PL Not Found !!! Please use correct VCF file!"<<std::endl;
            continue;
        }
        std::vector<std::string> pls = split(variant.FORMATS["PL"][s], ',');
        // if (pls.size() != (alts.size() + 1) * (alts.size() + 2) / 2)
        // {
        //     std::cerr<<"PL not Found"<<std::endl;
        //     exit(0);
        // }
    
        /// GATK: Even if GT is OK, PL still could be ".", so allele -> '?'
        if (pls.size() <3)
            continue;
        std::vector<float> PLs;
        // vectorized string to float
        std::transform(pls.begin(), pls.end(), std::back_inserter(PLs),
                        [](std::string &s) { return std::stof(s); });

        float minScore(1000000.0);
        for ( unsigned int p = 0; p < pls.size(); p++)
        {
            // GT's PL is PLs[ind]
            if (p == ind)
                continue;
            // float _het = PLs[p] - PLs[ind];
            minScore = std::min(minScore, PLs[p] - PLs[ind]);
        }
        // PL = -10 * \log{P(Genotype | Data)}, so that the PL value of the most likely genotype is 0,
        // PL_other - PL_gt >= heteroThereshold
        // PL: the lower, the more reliable            
        if (minScore >= opts->phredLikelihoodDifference)
        {
            if (GTs[0] > 0) 
            {
                if ((alts[GTs[0] - 1]) != "*") // GATK has *, covert to '?'
                {
                    std::string c = alts[GTs[0] - 1]; // do a copy
                    alleles[s] = c[0]; // string to char
                }
                hasAlt[GTs[0] - 1] = 1;
                //sampleAlts[s] = 1;
            } else {
                alleles[s] = (char) variant.REF[0];
            }
        }
    }
    return true;
}
bool VCF::parseStructralVariant(std::string& alleles, std::vector<std::string>& alts, std::vector<int>& hasAlt) 
{
    
    for (int s=0; s < strains.size(); s++)
    {
        std::string gt = variant.FORMATS["GT"][s];
        if (gt == "./." || gt == ".|.")
            continue;  // allele == '?'
        // replace
        std::replace(gt.begin(), gt.end(), '|', '/');
        std::vector<std::string> gts = split(gt, '/');
        
        if (gts.size() != 2)
            continue;
        if (gts[0] != gts[1])
            continue; // allele == '?'

        // genotype 0/0, 0/1, 1/1, 2/2 
        std::vector<int> GTs;
        // string to int, vectorize
        std::transform(gts.begin(), gts.end(), std::back_inserter(GTs),
                        [](std::string &s) { return std::stoi(s); });
        //unsigned int ind = (GTs[0] + 1) * (GTs[0] + 2) / 2 - 1;          

        if (GTs[0] > 0) 
        {
            if ((alts[GTs[0] - 1]) != "*") 
                alleles[s] = 'G';
            hasAlt[GTs[0] - 1] = 1;
        } else {
            alleles[s] = 'A'; // ref
        }
        
    }
    return true;
} 


void VCF::writeSNP(std::string& alleles)
{
    output <<"SNP_"<<variant.CHROM<<"_"<<variant.POS;
    output <<"\t"<<variant.CHROM<<"\t"<<variant.POS<<"\t";
    output <<variant.REF; // write REF

    // write allele pattern
    if (!samples.empty())
    {
        for(auto &s: samples)
        {
            int strIdx = strains.indexOf(s);
            output<<alleles[strIdx];
        }
    } else 
    {
        // for (auto &a: alleles) 
        output<<alleles;
    }
    output<<std::endl;
}

void VCF::writeTPED(std::string &alleles)
{
    std::string vartype;
    switch (variant.TYPE) {
        case 1:
            vartype = "SNP";
            break;
        case 2:
            vartype = "INDEL";
            break;
        case 3:
            vartype = "SV";
            break;
        default:
            vartype = "UNKNOWN";
    }
    tped <<variant.CHROM<<" "<<vartype<<"_"<<variant.CHROM<<"_"<<variant.POS;
    tped <<" 0 "<<variant.POS<<" ";
    if (variant.TYPE == 1) tped <<variant.REF<<" "<<variant.REF; // write REF
    if (variant.TYPE > 1)  tped <<"A A"; // convert to A as ref

    // write allele pattern
    if (!samples.empty())
    {
        for(auto &s: samples)
        {
            int strIdx = strains.indexOf(s);
            char a = alleles[strIdx];
            if (a == '?') a = '0';
            tped<<" "<<a<<" "<<a;
        }
    } else 
    {
        for (auto a: alleles) 
        {
            if (a=='?') a = '0';
            tped<<" "<<a<<" "<<a;
        }
    }
    tped<<std::endl;
}

void VCF::writeTFAM()
{
    // write allele pattern
    if (!samples.empty())
    {
        tfam<<"0 "<<outHeader<<" 0 0 1 0"<<std::endl;
        for (size_t i = 0; i < samples.size(); i ++)
        {
            tfam<<"0 "<<samples[i]<<" 0 0 1 "<<(i+1)<<std::endl;
        }
    } else 
    {
        tfam<<"0 "<<outHeader<<" 0 0 1 0"<<std::endl;
        for (int i=0; i < strains.size(); i++) 
            tfam<<"0 "<<strains.eltOf(i)<<" 0 0 1 "<<(i+1)<<std::endl;
    }
}


void VCF::writeStructralVariant(std::string &alleles)
{   
    std::string vartype;
    switch (variant.TYPE) {
        case 1:
            vartype = "SNP";
            break;
        case 2:
            vartype = "INDEL";
            break;
        case 3:
            vartype = "SV";
            break;
        default:
            vartype = "UNKNOWN";
    }
    std::unordered_map<std::string, std::string> INFO = variant.getINFO();
    output << vartype <<"_"<<variant.CHROM<<"_"<<variant.POS;
    if (variant.TYPE == 3) output <<"_"<<INFO["END"]<<"_"<<INFO["SVTYPE"];                
    output <<"\t"<<variant.CHROM<<"\t"<<variant.POS<<"\t";
    output << "A"; // write REF
    // write allele pattern
    if (!samples.empty())
    {
        for(auto &s: samples)
        {
            int strIdx = strains.indexOf(s);
            output<<alleles[strIdx];
        }
    } else 
    {
        // for (auto &a: alleles) 
        output<<alleles;
    }
    output<<std::endl;
            
}

void VCF::parseRecords() 
{
    int lineNum = 0;
    while (std::getline(*this->input, this->line))  
    {
        lineNum ++;
        // parse record
        variant = Variant(line);

        if (variant.QUAL < opts->qual) continue;
        if (CHROMOSOMES.find(variant.CHROM) == CHROMOSOMES.end()) 
        {
            std::cerr<<"chrom "<< variant.CHROM<< "is not recognized. Skip line: # "<<lineNum<<std::endl;
            continue;
        }

        std::string alleles(strains.size(), '?'); // allele pattern        
        std::vector<std::string> alts = split(variant.ALT, ',');
        std::vector<int> hasAlt(alts.size(), 0); // number of alternates
        std::vector<int> sampleAlts(strains.size(), 0); // sample is alt or not

        if (variant.TYPE == 1 && alts.size() > 1) // snp only allow 1 alternate
            continue;
        // string find not found, skip
        // std::unordered_map<std::string, std::string> INFO = variant.getINFO();

        if (variant.TYPE == 1 && std::strcmp(opts->variantType, "SNV") == 0)
        {
            bool ret = parseSNP(alleles, alts, hasAlt);
            if (!ret) continue;
            int numGoodAlt = std::accumulate(hasAlt.begin(), hasAlt.end(), 0);
            if (numGoodAlt == 1)
            {
                writeSNP(alleles);
                if (opts->plink) writeTPED(alleles);
            }
        } 
        else if ( variant.TYPE > 1 && (std::strcmp(opts->variantType, "SV") == 0 || std::strcmp(opts->variantType, "INDEL") == 0 ))
        {

            bool ret = parseStructralVariant(alleles, alts, hasAlt);
            if (!ret) continue;
            int numGoodAlt = std::accumulate(hasAlt.begin(), hasAlt.end(), 0);
            if (numGoodAlt == 1)
            {
                writeStructralVariant(alleles);
                if (opts->plink) writeTPED(alleles);
            }
        } 
    }
}