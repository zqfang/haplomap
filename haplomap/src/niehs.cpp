//
// Created by Zhuoqing Fang on 11/8/20.
//
#include <getopt.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include "utils.h"
#include "haplomap.h"

struct NIEHSOptions
{
    char *inputFileName;
    char *outputFileName;
    float heteroThreshold;
    float qual;
    bool verbose;

    // constructor
    NIEHSOptions() : inputFileName(nullptr), outputFileName(nullptr),
    heteroThreshold(20.0), qual(50.0), verbose(false){};
};



std::shared_ptr<NIEHSOptions> parseNIEHSOptions(int argc, char **argv) {
    int c;

    std::shared_ptr<NIEHSOptions> opts = std::make_shared<NIEHSOptions>();
    static struct option long_options_niehs[] = {
            {"help",          no_argument, nullptr,                           'h'},
            {"verbose",       no_argument, nullptr,                           'v'},
            {"input",         required_argument, nullptr,                           'i'},
            {"output",        required_argument, nullptr,                           'o'},
            {"hetero_thresh", optional_argument, nullptr, 't'},
            {"qual",          optional_argument, nullptr, 'q'},
            {nullptr, no_argument, nullptr,                           0}};

    const char *usage = "Convert VCF to NIEHS compact format)\n"
                        "\nusage: vcf2niehs [options]\n"
                        "\nrequired arguments:\n"
                        "    -i, --input           input VCF file, or 'stdin'\n"
                        "    -o, --output          output NIEHS compact file\n"
                        "\noptional arguments:\n"
                        "    -t, --hetero_thresh   heterozygous threshold, default 20.\n"
                        "    -q, --qual            qual field of VCF file, default 50.\n"
                        "    -v, --verbose\n"
                        "    -h, --help\n";

    while (true)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "hvi:o:t:q:", long_options_niehs, &option_index);

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
            case 'i':
            {
                opts->inputFileName = optarg;
                break;
            }

            case 'o':
            {
                // cout << "option -o with arg " << optarg << endl;
                opts->outputFileName = optarg;
                break;
            }

            case 't':
            {
                if (optarg == nullptr)
                    opts->heteroThreshold = 50.0;
                else
                    opts->heteroThreshold = std::stof(optarg);
                break;
            }
            case 'q':
                if (optarg == nullptr)
                    opts->qual = 20.0;
                else
                    opts->qual = std::stof(optarg);
                break;

            case '?':
            {
                /* getopt_long already printed an error message. */
                break;
            }
            default:
                abort();
        }
    }
    if (argc == 1)
    {
        std::cout<<usage<<std::endl;
        exit(1);
    }

    // if (opts->inputFileName == nullptr)
    // {
    //     std::cout<<usage<<std::endl;
    //     cout << "Required arg missing: input file name (-i)" << endl;
    //     exit(1);
    // }
    if (opts->outputFileName == nullptr)
    {
        std::cout<<usage<<std::endl;
        cout << "Required arg missing: output file name (-o)" << endl;
        exit(1);
    }

    // Print any remaining command line arguments (not options).
    if (optind < argc)
    {
        cout << "Extraneous things on command line: ";
        while (optind < argc)
        {
            cout << argv[optind++] << endl;
        }
        exit(1);
    }
    return opts;
}


struct Variant {
    std::string chrom;
    int pos;
    std::string ID;
    std::string ref;
    std::string alt;
    float qual;
    std::string filter;
    std::unordered_map<std::string, std::string> infos;
    std::unordered_map<std::string, std::vector<std::string>> formats;
    Variant() = default;
    explicit Variant(const std::string & record)
    {
        std::vector<std::string> rec = split(trim(record," \n\r"), '\t');
        chrom = rec[0];
        pos = std::stoi(rec[1]);
        ID = rec[2];
        ref = rec[3];
        alt = rec[4];
        qual = std::stof(rec[5]);
        filter = rec[6];

        std::vector<std::string> info = split(rec[7],';');
        std::vector<std::string> _info;
        for (auto &item: info) {
            _info = split(item,'=');
            // add values to hashmap
            infos[_info[0]] = _info.size() > 1 ? _info[1]: "true";
        }

        // format dict
        std::vector<std::string> format = split(rec[8], ':');
        std::vector<std::string> samples;
        for (size_t i = 9; i < rec.size(); i ++) {
            samples = split(rec[i],':');
            for (size_t j=0; j < format.size(); j++)
                // add values to hashmap
                formats[format[j]].push_back(samples[j]);
        }
    }
};

// Main for nhaploblocks
int main_niehs(int argc, char **argv)
{
    std::shared_ptr<NIEHSOptions> opts = parseNIEHSOptions(argc, argv);
    std::string rawHeader;
    std::string outHeader = "C57BL/6J";
    // output file
    std::ofstream output;
    output.open(opts->outputFileName);
    // read input file
    std::string line;
    Variant variant;
    std::vector<std::string> strains;
    std::istream* input;
    // std::ifstream input(opts->inputFileName);
    // stdin or ifstream
    if (opts->inputFileName == nullptr)
        input = &std::cin;  
    else
        input = new std::ifstream(opts->inputFileName, ios::in); 
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
        exit (1);
    }
    input->ignore(); // for clearing newline in cin

    std::cout<<"Parsing VCF file"<<std::endl;
    while (std::getline(*input, line))
    {
        // starts with "#"
        if (line.find('#') == 0) {
            if (line.find("#CHROM") == 0) {
                std::vector<std::string> header = split(trim(line," \n"),'\t');
                strains = std::vector<std::string>(header.begin()+9, header.end());
                output<<outHeader;
                for (auto &s: strains)
                    output <<"\t"<<s;
                output<<std::endl;
            }
            continue;
        }
        // parse record
        variant = Variant(line);

        if (variant.qual < opts->qual)
            continue;
        // skip indels, check if key is present
        if (variant.ref.size() > 1 || variant.infos.find("INDEL") != variant.infos.end())
            continue;

        //std::string alleles("?", strains.size());
        std::vector<std::string> alleles(strains.size(),"?");
        std::vector<std::string> alts = split(variant.alt, ',');
        std::vector<int> hasAlt(alts.size(), 0);

        for (size_t s=0; s < strains.size(); s++)
        {
            std::string gt = variant.formats["GT"][s];
            if (gt == "./." || gt == ".|.")
                continue;  // allele == '?'
            // replace
            std::replace(gt.begin(), gt.end(), '|', '/');
            std::vector<std::string> gts = split(gt, '/');

            if (gts.size() != 2)
                continue;
            if (gts[0] == "." || gts[0] != gts[1])
                continue; // allele == '?'

            // genotype 0/0, 0/1, 1/1
            std::vector<int> GTs;
            // string to int, vectorize
            std::transform(gts.begin(), gts.end(), std::back_inserter(GTs),
                            [](std::string &s) { return std::stoi(s); });
            unsigned int ind = (GTs[0] + 1) * (GTs[0] + 2) / 2 - 1;

            if (variant.formats.find("PL") == variant.formats.end())
            {
                std::cerr<<"PL Not Found !!! Please use correct VCF file!"<<std::endl;
                continue;
            }
            std::vector<std::string> pls = split(variant.formats["PL"][s], ',');
//                if (pls.size() != (alts.size() + 1) * (alts.size() + 2) / 2)
//                {
//                    std::cerr<<"PL not Found"<<std::endl;
//                    exit(0);
//                }
            float minScore(1000000000.0);
            std::vector<float> PLs;
            std::transform(pls.begin(), pls.end(), std::back_inserter(PLs),
                            [](std::string &s) { return std::stof(s); });
            // homozygous - heterozyous > heteroThereshold
            for ( unsigned int p = 0; p < PLs.size(); p++)
            {
                if (p == ind)
                    continue;
                if ((PLs[p] - PLs[ind]) < minScore)
                    minScore = PLs[p] - PLs[ind];
            }

            if (minScore >= opts->heteroThreshold)
            {
                if (GTs[0] > 0) {
                    alleles[s] = alts[GTs[0] - 1];
                    hasAlt[GTs[0] - 1] = 1;
                } else {
                    alleles[s] = variant.ref;
                }
            }
        }

        // find good alt
        int numGoodAlt = std::accumulate(hasAlt.begin(), hasAlt.end(), 0);
        if (numGoodAlt == 1)
        {
            output <<"SNP_"<<variant.chrom<<"_"<<variant.pos;
            output <<"\t"<<variant.chrom<<"\t"<<variant.pos<<"\t";
            output <<variant.ref;
            for (auto &a: alleles)
                output<<a;
            output<<std::endl;
        }

    }
    
    // input.close();
    output.close();
    if (opts->inputFileName != nullptr)
        delete input;
    std::cout<<"Job done."<<std::endl;
    return 0;
}