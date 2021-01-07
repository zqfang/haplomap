//
// Created by Zhuoqing Fang on 11/8/20.
//

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <getopt.h>
// #include <zlib.h>
#include "utils.h"
#include "haplomap.h"
#include "dynum.h"

#define GZ_BUF_SIZE 1048576

struct NIEHSOptions
{
    char *inputFileName;
    char *outputFileName;
    char *sampleFileName;
    float phredLikelihoodDifference; // PL_other - PL_homo > 20
    float qual;  // QUAL field
    
    /// allele depth of coverage (AD)
    /// variants < AD set to be LowQual
    float alleleDepth;
    /// see here https://samtools.github.io/bcftools/howtos/variant-calling.html
    /// ratio = %MAX(AD) / %MAX(DP). variants < ratio set to be LowQual
    float ratio; 
    bool verbose;

    // constructor
    NIEHSOptions() : inputFileName(nullptr), outputFileName(nullptr), 
                     sampleFileName(nullptr), phredLikelihoodDifference(20.0), 
                     qual(50.0), alleleDepth(3.0), ratio(0.1), verbose(false) {};
};



std::shared_ptr<NIEHSOptions> parseNIEHSOptions(int argc, char **argv) 
{
    std::shared_ptr<NIEHSOptions> opts = std::make_shared<NIEHSOptions>();
    static struct option long_options_niehs[] = {
            {"help",           no_argument, nullptr,       'h'},
            {"verbose",        no_argument, nullptr,       'v'},
            {"input",          required_argument, nullptr, 'i'},
            {"output",         required_argument, nullptr, 'o'},
            {"samples",        optional_argument, nullptr, 's'},
            {"pl-diff",        optional_argument, nullptr, 'p'},
            {"qual",           optional_argument, nullptr, 'q'},
            {"allele-depth",   optional_argument, nullptr, 'a'},
            {nullptr,          no_argument, nullptr,        0}};

    const char *usage = "Convert VCF to NIEHS compact format\n"
                        "\nusage: vcf2niehs [options] <in.vcf> \n"
                        "\nrequired arguments:\n"
                        "    in.vcf                Input sorted VCF file or stdin\n"
                        "    -o, --output          Output file name\n"
                        "\noptional arguments:\n"
                        "    -s, --samples         New sample order. One name per line.\n"
                        "    -p, --pl-diff         Phred-scaled genotype likelihood (PL) difference. Default 20.\n"
                        "                          GT's PL must at least pl-diff unit lower than any other PL value. \n"
                        "                          The larger, the more confident \n"
                        "    -q, --qual            QUAL field of VCF file. Only keep variant > qual. Default 50. \n"
                        "    -a, --allele-depth    Min allele depth (AD) of samples. Default 3.\n"
                        "    -r, --ratio           Min ratio (%MAX(AD) / %MAX(DP)) of samples. Default 0.1.\n"
                        "    -v, --verbose\n"
                        "    -h, --help\n";
    if (argc == 1)
    {
        std::cout<<usage<<std::endl;
        exit(1);
    }
    
    
    int c;
    while (true)
    {

       int option_index = 0;
       c = getopt_long(argc, argv, "hv:a:o:p:q:r:s:p:", long_options_niehs, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
        {
            break;
        }

        switch (c)
        {

            case 'h':
            {
                cout << usage << endl;
                exit(0);
                //break;
            }

            case 'v':
            {
                opts->verbose = true;
                break;
            }

            case 'o':
            {
                // cout << "option -o with arg " << optarg << endl;
                opts->outputFileName = optarg;
                break;
            }
            case 's':
            {
                opts->sampleFileName = optarg;
                break;
            }

            case 'p':
            {
                if (optarg != nullptr)
                    opts->phredLikelihoodDifference = std::stof(optarg);
                break;
            }
            case 'q':
            {
                if (optarg != nullptr)
                    opts->qual = std::stof(optarg);
                break;
            }
            case 'a':
            {
                if (optarg != nullptr)
                    opts->alleleDepth = std::stof(optarg);
                break;
            }  
            case 'r':
            {
                if (optarg != nullptr)
                    opts->ratio = std::stof(optarg);
                break;
            }           
            case '?':
            {
                /* getopt_long already printed an error message. */
                // this is for unkwon options
                break;
            }
            default:
                abort();
        }
    }

    if (opts->outputFileName == nullptr)
    {
        std::cout<<usage<<std::endl;
        std::cout << "Required argument missing: --output (-o)" << std::endl;
        exit(1);
    }
    // Print any remaining command line arguments (not options).
    if (optind < argc)
    {
        while (optind < argc)
        {
            opts->inputFileName = argv[optind++];
        }
    }
    return opts;
}


class Variant {

public:
    bool isINDEL;
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
    Variant(std::string record)
    {
        std::vector<std::string> rec = split(trim(record), '\t');
        CHROM = rec[0];
        POS = std::stoi(rec[1]);
        ID = rec[2];
        REF = rec[3];
        ALT = rec[4];
        QUAL = std::stof(rec[5]);
        FILTER = rec[6];
        INFO = rec[7];
        isINDEL = isIndel();
        // format dict
        parseFORMAT(rec);
       
    }
    std::unordered_map<std::string, std::string> getINFO() 
    {
        std::unordered_map<std::string, std::string> infos;
        std::vector<std::string> info = split(INFO,';');
        std::vector<std::string> _info;
        for (auto &item: info) {
            _info = split(item,'=');
            // add values to hashmap
            infos[_info[0]] = _info.size() == 2 ? _info[1]: "true";
        }    
        return infos;
    }
    float maxAlleleDepth()
    {
        if (this->FORMATS.find("AD") == this->FORMATS.end())
            return 0; // if AD not found, then we should skip 
        
        float AD(0);
        for (auto & a: this->FORMATS["AD"])
        {
            // AD is a vector for gt (ref0, alt1, alt2)
            std::vector<std::string> ad = split(a, ',');
            for (auto it = ad.begin() + 1; it != ad.end(); it++) 
            {
                if (*it != ".")
                    AD = std::max(std::stof(*it), AD);
            }
        }           
        return AD;   
    }

    float maxDepth()
    {
        if (this->FORMATS.find("DP") == this->FORMATS.end() )
            return 0; // if DP not found, then we should skip 
        
        float DP(0);
        /// FIXME: if DP is .
        for (auto & a: this->FORMATS["DP"])
        {   
            if (a != ".")
                DP = std::max(std::stof(a), DP);
        }
        return DP;
    }
        
private:
    void parseFORMAT(const std::vector<std::string> & rec) {
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
    bool isIndel()
    {
        if (this->REF.size() > 1)
            return true;
        std::vector<std::string> alts = split(this->ALT, ',');
        for (auto &alt: alts)
        {
            if (alt.size() > 1)
                return true;
        }
        return false;
    }
    
    
};


// gzip file read

// bool gzLoad(char* gzfn, std::string &out)
// {
// 	//open .gz file
// 	gzFile gzfp = gzopen(gzfn,"rb");
// 	if(!gzfp)
// 	{
// 		return false;
// 	}

// 	//read and add it to out
// 	unsigned char buf[GZ_BUF_SIZE];
// 	int have;
// 	while( (have = gzread(gzfp,buf,GZ_BUF_SIZE)) > 0)
// 	{
// 		out.append((const char*)buf,have);
// 	}

// 	//close .gz file
// 	gzclose(gzfp);
// 	return true;
// }


// Parsing VCF
int main_niehs(int argc, char **argv)
{
    std::shared_ptr<NIEHSOptions> opts = parseNIEHSOptions(argc, argv);
    // first allele is the reference
    std::string outHeader = "C57BL/6J";
    std::string line;
    Variant variant;
    Dynum<std::string> strains;
    std::vector<std::string> samples;

    // read input file
    std::istream* input;
    // std::ifstream input(opts->inputFileName);
    // stdin or ifstream
    if (opts->inputFileName == nullptr)
    {
        std::cout<<"Read VCF: from stdin"<<std::endl;
        input = &std::cin;
    }
    else
    {
        std::cout<<"Read VCF: "<< opts->inputFileName << std::endl;
        input = new std::ifstream(opts->inputFileName, ios::in); 
    }
    // check status
    if ( input->fail() ) 
    {
        std::cerr << "Error: The requested file (" 
                << opts->inputFileName
                << ") " 
                << "could not be opened. "
                << "Error message: ("
                << std::strerror(errno)
                << "). Exiting!" << std::endl;
        exit(1);
    }
    input->ignore(); // for clearing newline in cin

    // output file
    std::ofstream output;
    output.open(opts->outputFileName);      
    // output header
    output<<outHeader;
    
    // read input sample names
    if (opts->sampleFileName != nullptr) 
    {
        std::cout << "Read Sample Names: "<< opts->sampleFileName << std::endl;
        std::ifstream sinput(opts->sampleFileName);
        while (std::getline(sinput, line))
        {
            // starts with "#"
            if (line.find('#') == 0) continue;
            samples.push_back(line);    
            output <<"\t"<<line;            
        } 
        output<<std::endl;
        sinput.close();
    } 
    

    // read vcf header
    std::cout<<"Parser Header"<<std::endl;
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
             
            if (opts->sampleFileName != nullptr) 
            {
                if (strains.size() < samples.size())
                {
                    std::cerr<<"Input (--samples) size larger than vcf sample size !!!"
                             <<std::endl;
                    exit(-1);
                }
                break;
            }
            // write header
            for (auto itr = header.begin()+9; itr != header.end(); itr++)
                output <<"\t"<< *itr;
            output<<std::endl;
            break;
        }  
    }
    // check samples in strains
    if (!samples.empty()) 
    {
        for (auto &s: samples)
        {
            int inStrain = strains.hasIndex(s);
            if (inStrain == -1)
            {
                std::cerr<<"Input name: "<<s<<"  NOT in VCF samples"<<std::endl;
                exit(1);
            }
        }
    }
    int lineNum = 0;
    // Now read all records
    std::cout<<"Parsing Variants"<<std::endl;
    while (std::getline(*input, line))
    {
        lineNum ++;
        // parse record
        variant = Variant(line);

        if (variant.QUAL < opts->qual)
            continue;
        float AD = variant.maxAlleleDepth();
        if (AD < opts->alleleDepth || (AD / variant.maxDepth()) < opts->ratio)
            continue; 
        // string find not found, skip
        if (variant.isINDEL || variant.INFO.find("INDEL") != std::string::npos) 
            continue;

        //std::string alleles("?", strains.size()); // not easy to do inplace replacement
        
        std::vector<std::string> alts = split(variant.ALT, ',');
        std::vector<int> hasAlt(alts.size(), 0); // number of alternates
        std::vector<int> sampleAlts(strains.size(), 0); // sample is alt or not

        // if (alts.size() > 1)
        //     continue;

        std::vector<std::string> alleles(strains.size(),"?");        
        for (size_t s=0; s < strains.size(); s++)
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
            // PL_other - PL_gt >= heteroThereshold
            // PL: the lower, the more reliable            
            if (minScore >= opts->phredLikelihoodDifference)
            {
                if (GTs[0] > 0) 
                {
                    if ((alts[GTs[0] - 1]) != "*") // GATK has *, covert to '?'
                        alleles[s] = alts[GTs[0] - 1];
                    hasAlt[GTs[0] - 1] = 1;
                    sampleAlts[s] = 1;
                } else {
                    alleles[s] = variant.REF;
                }
            }
        }
        // find good alt
        int numGoodAlt = std::accumulate(hasAlt.begin(), hasAlt.end(), 0);
        int sampleGoodAlt = std::accumulate(sampleAlts.begin(), sampleAlts.end(), 0);
        /// FIXME: if there are 2 alts, and all input samples are not ref, do we still keep it? yes
        if (numGoodAlt == 1 || (sampleGoodAlt == sampleAlts.size() && numGoodAlt == 2)) /// FIXED: 
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
                for (auto &a: alleles) 
                    output<<a;
            }
            output<<std::endl;
        }
        // if (lineNum % 100000 == 0)
        //     std::cout<<"Parsed line #: "<<lineNum<<std::endl;

    }
    
    // input.close();
    output.close();
    if (opts->inputFileName != nullptr)
        delete input;
    std::cout<<"Job done."<<std::endl;
    return 0;
}
