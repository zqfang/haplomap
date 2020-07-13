//
// Created by Zhuoqing Fang on 7/12/20.
//
#include <memory>
#include <iostream>
#include <fstream>
#include "gsl/gsl_matrix.h"
#include "haplomap.h"
#include "eigen.h"


struct EigenOptions
{
    unsigned L;
    char *inputFileName;
    char *outputFileName;
    bool verbose;

    // constructor
    EigenOptions() : L(4), inputFileName(NULL), outputFileName(NULL), verbose(false){};
};

std::shared_ptr<EigenOptions> parseEigenOptions(int argc, char **argv)
{
    int c;

    std::shared_ptr<EigenOptions> opts(new EigenOptions());

    static struct option long_options_ghmap[] = {
            {"help", no_argument,            0,                          'h'},
            {"verbose", no_argument,         0,                          'v'},
            {"input", required_argument,     0,                          'i'},
            {"output", required_argument,    0,                          'o'},
            {"dimension", optional_argument, reinterpret_cast<int *>(4), 'd'},
            {0, 0,                           0,                          0}};

    const char *usage = "Perform reduction on the data dimension (rows)\n"
                        "\nusage: pca [options]\n"
                        "\nrequired arguments:\n"
                        "    -i, --input           input file (M x N matrix)\n"
                        "    -o, --output          output file (L x N matrix)\n"
                        "\noptional arguments:\n"
                        "    -d, --dimension       dimensions of reduction, default 4.\n"
                        "    -v, --verbose\n"
                        "    -h, --help\n";

    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "hvi:o:d:", long_options_ghmap, &option_index);

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
                break;
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

            case 'd':
            {
                // cout << "option -p with arg " << optarg << endl;
                opts->L = std::atoi(optarg);
                break;
            }

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

    if (NULL == opts->inputFileName)
    {
        std::cout<<usage<<std::endl;
        cout << "Required arg missing: input file name (-i)" << endl;
        exit(1);
    }
    else if (NULL == opts->outputFileName)
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


int main_eigen(int argc, char **argv) {
    std::shared_ptr<EigenOptions> opts = parseEigenOptions(argc, argv);
    std::shared_ptr<EigenMat> eigens;
    eigens = std::make_shared<EigenMat>(opts->inputFileName, false, "\t");
    gsl_matrix *results = eigens->pca(opts->L);

    // write output file
    std::ofstream output;
    output.open(opts->outputFileName);
    if (output.is_open()) {
        for (size_t i = 0; i < results->size1; ++i) {
            for (size_t j = 0; j < results->size2; ++j)
                output << gsl_matrix_get(results, i, j) << "\t";
            output << std::endl;
        }
    }
    output.close();
    gsl_matrix_free(results);
    std::cout<<"Job done."<<std::endl;
    return 0;
}
