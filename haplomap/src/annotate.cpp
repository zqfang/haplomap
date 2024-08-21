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
    bool maxPrioritize;
    bool verbose;

    VEPOptions(): inputVEPName(nullptr), inputStrainName(nullptr), 
                    outputSUMName(nullptr), outputCSQName(nullptr),
                    variantType((char*)"all"), maxPrioritize(false), 
                    verbose(false) {}
};

std::shared_ptr<VEPOptions> parseVEPOptions(int argc, char **argv) 
{
    std::shared_ptr<VEPOptions> opts = std::make_shared<VEPOptions>();
    // see summary section of how to set the parameter: https://azrael.digipen.edu/~mmead/www/Courses/CS180/getopt.html
    static struct option long_options_vep[] = {
            {"help",           no_argument, nullptr,       'h'},
            {"verbose",        no_argument, nullptr,       'v'},
            //{"input",          optional_argument, nullptr, 'i'}, // 
            {"output",         required_argument, nullptr, 'o'},
            {"csq",            required_argument, nullptr, 'c'},
            {"samples",        required_argument, nullptr, 's'},
            {"type",           required_argument, nullptr, 't'},
            {"prioritize",     required_argument, nullptr, 'p'},
            {nullptr,          no_argument, nullptr,        0}};

    const char *usage = "Convert ensembl-vep to eblocks (-g) input\n"
                        "\nUsage: annotate [options] <in.vep.txt> \n"
                        "\nRequired arguments:\n"
                        "    in.vep.txt             Input ensembl-VEP tab format file name\n"
                        "    -o, --output           Output file name, for (eblocks -g)\n"
                        "\nOptional arguments:\n"
                        "    -c,  --csq             Output a annotation with impact score and sample names.\n"
                        "    -s,  --samples         Only write annotation for the input samples (e.g. eblocks -s).\n"
                        "    -t,  --type            Select variant type: [snp|indel|sv|all]. Default: all\n"
                        "    -p,  --prioritize      Whether aggregate variant annotation by max impact score. Default: false"
                        "    -v,  --verbose\n"
                        "    -h,  --help\n"
                        "\nImportant message:\n"
                        "Please run ensemble-vep containing following flags: \n"
                        "    vep --fasta --individual_zyg all --everything\n\n"
                        "Structural variant input format for ensembl-vep, please ref to: \n"
                        "    https://ensembl.org/info/docs/tools/vep/vep_formats.html#sv \n";

    if (argc == 1)
    {
        std::cout<<usage<<std::endl;
        exit(1);
    }
    
    
    int c;
    while (true)
    {

       int option_index = 0;
       // : -> expect an associated value
       c = getopt_long(argc, argv, "hvpc:s:o:t:", long_options_vep, &option_index);

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
            case 'p':
            {
                opts->maxPrioritize = true;
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
    if (opts->verbose) 
    { 
        std::cout<<"Read ensembl-vep results ..."<<std::endl;
        // std::count<<"vep --fasta --individual_zyg all --everything"<<std::endl;
    }
    
    // handle input variant type
    // to lower case
    std::string varType(opts->variantType);
    std::transform(varType.begin(), varType.end(), varType.begin(), ::tolower);
    // C-style
    // char* varType = strdup(opts->variantType);
    // char* name = varType;
    // while (*name) 
    // {
    //     *name = std::tolower(*name);
    //     name++;
    // }
    // ...
    // free(varType)
    std::vector<std::string> v = {"snp","snv", "indel", "sv", "all"};
    if (std::find(v.begin(), v.end(), varType) == v.end())
    {
        std::cerr<<"Variant type (-t) error. Input one of these: snp, indel, sv"<<std::endl;
        std::exit(1);
    }
    if (varType == "snp") varType = "snv";
    opts->variantType = (char*)varType.c_str();
    // read data
    vep.readVEP(opts->inputVEPName, (char*)"\t", (char*)varType.c_str());
    // write
    if (opts->verbose) std::cout<<"Write Variant Annotation for eblocks"<<std::endl;
    if (opts->verbose && opts->maxPrioritize) 
    {
        std::cout<<"    Write Annotation with max prioritization procedure"<<std::endl;
        std::cout<<"    For each variant (row), only write the most impactful variant per gene.\n";
        std::cout<<"    If multiple variant consequence has the same impact score, select one randomly.\n";
    }
    vep.writeVEPCsq(opts->outputCSQName, opts->maxPrioritize);
    if (opts->outputSUMName != nullptr)
    {
        if (opts->verbose) std::cout<<"Write CodonFlag Annotation for eblocks"<<std::endl;
        vep.writeVEPImpact(opts->outputSUMName);
    }
    if (opts->verbose) std::cout<<"Job done"<<std::endl;
    return 0;
}