#include <memory>
#include <iostream>
#include <string>
#include <getopt.h>
#include "haplomap.h"
#include "vep.h"


struct VEPOptions {
    char *inputVEPName;
    char *inputStrainName;
    char *outputSUMName;
    char *outputCSQName;
    char *variantType;
    bool verbose;

    VEPOptions(): inputVEPName(nullptr), inputStrainName(nullptr), 
                    outputSUMName(nullptr), outputCSQName(nullptr),
                    variantType((char*)"all"), verbose(false) {}
};

std::shared_ptr<VEPOptions> parseVEPOptions(int argc, char **argv) 
{
    std::shared_ptr<VEPOptions> opts = std::make_shared<VEPOptions>();
    static struct option long_options_niehs[] = {
            {"help",           no_argument, nullptr,       'h'},
            {"verbose",        no_argument, nullptr,       'v'},
            {"input",          required_argument, nullptr, 'i'},
            {"output",         required_argument, nullptr, 'o'},
            {"csq",            optional_argument, nullptr, 'c'},
            {"samples",        optional_argument, nullptr, 's'},
            {"type",           optional_argument, nullptr, 't'},
            {nullptr,          no_argument, nullptr,        0}};

    const char *usage = "Convert ensembl-VEP to eblocks (-g) input\n"
                        "\nUsage: annotate [options] <in.vcf> \n"
                        "\nRequired arguments:\n"
                        "    in.vep.txt             Input sorted ensembl-VEP file or stdin\n"
                        "    -o, --output           Output file name, for (eblocks -g)\n"
                        "\nOptional arguments:\n"
                        "    -c,  --csq             Output a VEP consequence annotation file.\n"
                        "    -s,  --samples         Only get annotation from these samples, same to (eblocks -s).\n"
                        "    -t,  --type            Select variant type: [snp|indel|sv|all]. Default: all\n"
                        "    -v,  --verbose\n"
                        "    -h,  --help\n"
                        "\nImportant message:\n"
                        "Please run ensemble-vep containing following flags: \n"
                        "    vep --fasta --format vcf --individual all --gencode_basic --everything\n";

    if (argc == 1)
    {
        std::cout<<usage<<std::endl;
        exit(1);
    }
    
    
    int c;
    while (true)
    {

       int option_index = 0;
       c = getopt_long(argc, argv, "hvc:s:o:t:", long_options_niehs, &option_index);

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

            case 'c':
            {
                // cout << "option -c with arg " << optarg << endl;
                if (optarg != nullptr)
                    opts->outputSUMName = optarg;
                break;
            }
            case 's':
            {
                opts->inputStrainName = optarg;
                break;
            }
            case 't':
            {
                if (optarg != nullptr)
                    opts->variantType = optarg;
                break;
            }  
            case 'o':
            {
                if (optarg != nullptr)
                    opts->outputCSQName = optarg;
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

    if (opts->outputCSQName == nullptr)
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
            opts->inputVEPName = argv[optind++];
        }
    }
    return opts;
}


int main_annot(int argc, char **argv)
{
    std::shared_ptr<VEPOptions> opts = parseVEPOptions(argc, argv);

    VarirantEeffectPredictor vep(opts->inputVEPName, opts->inputStrainName);
    if (opts->verbose) std::cout<<"Read ensembl-vep Annotation"<<std::endl;
    
    // handle input variant type
    // to lower case
    char* name = opts->variantType;
    while (*name) 
    {
        *name = tolower(*name);
        name++;
    }
    std::vector<std::string> v = {"snp","snv", "indel", "sv", "all"};
    if ((opts->variantType != nullptr ) && 
        (std::find(v.begin(), v.end(), std::string(opts->variantType)) == v.end()))
    {
        std::cerr<<"Variant type (-t) error. Input one of these: snp, indel, sv"<<std::endl;
        std::exit(1);
    }

    if (std::strcmp(opts->variantType, "snp") == 0)
    {
        opts->variantType = (char*)"snv";
    }
    // read data
    vep.readVEP(opts->inputVEPName, (char*)"\t", opts->variantType);
    // write
    if (opts->verbose) std::cout<<"Write detail Annotation"<<std::endl;
    vep.writeVEPCsq(opts->outputCSQName);
    if (opts->outputSUMName != nullptr)
    {
        if (opts->verbose) std::cout<<"Write CodonFlag Annotation"<<std::endl;
        vep.writeVEPImpact(opts->outputSUMName);
    }
    if (opts->verbose) std::cout<<"Job done"<<std::endl;
    return 0;
}