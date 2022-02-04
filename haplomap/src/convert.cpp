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
                        "\nusage: convert [options] <in.vcf> \n"
                        "\nrequired arguments:\n"
                        "    in.vcf                Input sorted VCF file or stdin\n"
                        "    -o, --output          Output file name\n"
                        "\noptional arguments:\n"
                        "    -s,  --samples         New sample order. One name per line.\n"
                        "    -t,  --type            Select variant type: [snp|indel|sv]. Default: snp\n"
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


// Parsing VCF
int main_convert(int argc, char **argv)
{
    std::shared_ptr<VCFOptions> opts = parseNIEHSOptions(argc, argv);
    // 
    // to lower case
    char* name = opts->variantType;
    while (*name) 
    {
        *name = tolower(*name);
        name++;
    }
    if (std::strcmp(opts->variantType, "snp") == 0)
    {
        opts->variantType = (char*)"snv";
    }
    std::vector<std::string> v = {"snp","snv", "indel", "sv"};
    if ((opts->variantType != nullptr ) && 
        (std::find(v.begin(), v.end(), std::string(opts->variantType)) == v.end()))
    {
        std::cerr<<"Variant type error. Input one of these: snp, indel, sv"<<std::endl;
        std::exit(1);
    }
    ///
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
