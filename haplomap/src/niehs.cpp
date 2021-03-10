//
// Created by Zhuoqing Fang on 11/8/20.
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
#include <getopt.h>
#include "utils.h"
#include "haplomap.h"
#include "dynum.h"
#include "vcf.h"



std::shared_ptr<VCFOptions> parseNIEHSOptions(int argc, char **argv) 
{
    std::shared_ptr<VCFOptions> opts = std::make_shared<VCFOptions>();
    static struct option long_options_niehs[] = {
            {"help",           no_argument, nullptr,       'h'},
            {"verbose",        no_argument, nullptr,       'v'},
            {"input",          required_argument, nullptr, 'i'},
            {"output",         required_argument, nullptr, 'o'},
            {"samples",        optional_argument, nullptr, 's'},
            {"type",           optional_argument, nullptr, 't'},
            {"pl-diff",        optional_argument, nullptr, 'p'},
            {"qual",           optional_argument, nullptr, 'q'},
            {"allelic-depth",  optional_argument, nullptr, 'a'},
            {"mapping-quality",optional_argument, nullptr, 'M'},
            {"strand-bias",    optional_argument, nullptr, 'S'},
            {"ratio",          optional_argument, nullptr, 'r'},
            {nullptr,          no_argument, nullptr,        0}};

    const char *usage = "Convert VCF to NIEHS compact format\n"
                        "\nusage: vcf2niehs [options] <in.vcf> \n"
                        "\nrequired arguments:\n"
                        "    in.vcf                Input sorted VCF file or stdin\n"
                        "    -o, --output          Output file name\n"
                        "\noptional arguments:\n"
                        "    -s,  --samples         New sample order. One name per line.\n"
                        "    -t,  --type            Select variant type: snps,indels,bnd,sv.\n"
                        "    -q,  --qual            QUAL field of VCF file. Only keep variant > qual. Default 50. \n"
                        "\nSNP only arguments:\n"
                        "    -p,  --pl-diff         Phred-scaled genotype likelihood (PL) difference. Default 20.\n"
                        "                           GT's PL must at least pl-diff unit lower than any other PL value. \n"
                        "                           The greater, the more '?' in the output. \n" 
                        "    -a,  --allelic-depth   Min allelic depth (AD) of samples. Default 3.\n"
                        "    -m,  --mapping-quality Min average mapping quality. Default 20. \n"
                        "    -b,  --strand-bias     Max Phred-scaled pvalue for strand bias (the lower, the better). Default 50. \n"
                        "    -r,  --ratio           Min ratio of (%MAX(AD) / %MAX(DP)). Default 0.1.\n\n"
                        "    -v,  --verbose\n"
                        "    -h,  --help\n";
    if (argc == 1)
    {
        std::cout<<usage<<std::endl;
        exit(1);
    }
    
    
    int c;
    while (true)
    {

       int option_index = 0;
       c = getopt_long(argc, argv, "hva:m:o:p:q:r:s:t:b:", long_options_niehs, &option_index);

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
            case 't':
            {
                if (optarg != nullptr)
                    opts->variantType = optarg;
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
            case 'm':
            {
                if (optarg != nullptr)
                    opts->mappingQuality = std::stof(optarg);
                break;
            }
            case 'b':
            {
                if (optarg != nullptr)
                    opts->strandBiasPhredPval = std::stof(optarg);
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



void parseStructralVariant(Variant& variant, std::vector<std::string> & alleles, 
                           std::vector<std::string>& alts, std::vector<int>& hasAlt,
                           Dynum<std::string> strains) 
{
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
        //unsigned int ind = (GTs[0] + 1) * (GTs[0] + 2) / 2 - 1;          

        if (GTs[0] > 0) 
        {
            if ((alts[GTs[0] - 1]) != "*") 
                alleles[s] = '1';
            hasAlt[GTs[0] - 1] = 1;
        } else {
            alleles[s] = '0';
        }
        
    }
    return ;
} 


// Parsing VCF
int main_niehs2(int argc, char **argv)
{
    std::shared_ptr<VCFOptions> opts = parseNIEHSOptions(argc, argv);
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
        if (opts->verbose)
            std::cout<<"Read VCF: from stdin"<<std::endl;
        input = &std::cin;
    }
    else
    {
        if (opts->verbose)
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
        if (opts->verbose)
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
    
    if (opts->verbose)
        std::cout<<"Parsing Header"<<std::endl;
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
    if (opts->verbose)
        std::cout<<"Parsing Variants"<<std::endl;
    while (std::getline(*input, line))
    {
        lineNum ++;
        // parse record
        variant = Variant(line);

        if (variant.QUAL < opts->qual)
            continue;

        //std::string alleles("?", strains.size()); // not easy to do inplace replacement
        
        std::vector<std::string> alts = split(variant.ALT, ',');
        std::vector<int> hasAlt(alts.size(), 0); // number of alternates
        std::vector<int> sampleAlts(strains.size(), 0); // sample is alt or not

        // if (alts.size() > 1) // only allow 1 alternate
        //     continue;
 
        // string find not found, skip

        std::unordered_map<std::string, std::string> INFO = variant.getINFO();
        std::vector<std::string> alleles(strains.size(),"?");  
        if (variant.isSNP && std::strcmp(opts->variantType, "snps") == 0)
        {      
            if (variant.REF == "N")
                continue;
            
            if (INFO["MQ"].size() < 1 || std::stof(INFO["MQ"]) < opts->mappingQuality)
                continue;
        
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
                continue; // filter strand bias variant 

            float AD = variant.maxAllelicDepth();
            if (AD < opts->alleleDepth || (AD / variant.maxDepth()) < opts->ratio)
                continue;
            // start parsing samples
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
        } else {// strutral variant
            parseStructralVariant(variant,  alleles, alts,  hasAlt, strains);
        }
        // find good alt
        int numGoodAlt = std::accumulate(hasAlt.begin(), hasAlt.end(), 0);
        //int sampleGoodAlt = std::accumulate(sampleAlts.begin(), sampleAlts.end(), 0);
        /// FIXME: if there are 2 alts, and all input samples are not ref, do we still keep it? NO.
        // if (numGoodAlt == 1 || (sampleGoodAlt == sampleAlts.size() && numGoodAlt == 2)) /// FIXED: 
        if (numGoodAlt == 1) /// FIXED:
        {
            
            if (variant.isSNP && std::strcmp(opts->variantType, "snps") == 0)
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
            else if (!variant.isSNP && std::strcmp(opts->variantType, "sv") == 0) 
            {   
                output <<"SNP_"<<variant.CHROM<<"_"<<variant.POS; 
                output <<"_"<<INFO["END"]<<"_"<<INFO["SVTYPE"];                
                output <<"\t"<<variant.CHROM<<"\t"<<variant.POS<<"\t";
                output << "0"; // write REF
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
        }
        // if (lineNum % 100000 == 0)
        //     std::cout<<"Parsed line #: "<<lineNum<<std::endl;

    }
    
    output.close();
    if (opts->inputFileName != nullptr)
        delete input;
    if (opts->verbose)
        std::cout<<"Job done."<<std::endl;
    return 0;
}


// Parsing VCF
int main_niehs(int argc, char **argv)
{
    std::shared_ptr<VCFOptions> opts = parseNIEHSOptions(argc, argv);

    std::string line;
    //Dynum<std::string> strains;
    std::vector<std::string> samples;
     
    // read input file
    std::istream* input;
    // std::ifstream input(opts->inputFileName);
    // stdin or ifstream
    if (opts->inputFileName == nullptr)
    {
        if (opts->verbose)
            std::cout<<"Read VCF: from stdin"<<std::endl;
        input = &std::cin;
    }
    else
    {
        if (opts->verbose)
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

    // read input sample names
    if (opts->sampleFileName != nullptr) 
    {
        if (opts->verbose)
            std::cout << "Read Sample Names: "<< opts->sampleFileName << std::endl;
        std::ifstream sinput(opts->sampleFileName);
        while (std::getline(sinput, line))
        {
            // starts with "#"
            if (line.find('#') == 0) continue;
            samples.push_back(line);    
        } 
        sinput.close();
    } 
    
    if (opts->verbose)
        std::cout<<"Parsing Header"<<std::endl;
        
    // read vcf header
    VCF vcf(input, opts, samples);
    // Now read all records
    if (opts->verbose)
        std::cout<<"Parsing Variants"<<std::endl;
    vcf.parseRecords();


    if (opts->inputFileName != nullptr)
        delete input;
    if (opts->verbose)
        std::cout<<"Job done."<<std::endl;
    return 0;
}
