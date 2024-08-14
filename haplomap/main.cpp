//
// Created by Zhuoqing Fang on 7/11/20.
//

#ifndef __HAPLOMAP_VER__
#define __HAPLOMAP_VER__ "0.1.2"
#endif

#if defined(__clang__)
#define __COMPILER__ "clang++"
#elif defined(WIN32)
#define ___COMPILER__ "MSVC"
#else
#define __COMPILER__ "g++"
#endif

#include <iostream>
#include <cstring>
#include "haplomap.h"

using namespace std;

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help;
} cmd_t;

static cmd_t cmds [] = {
        {
            .func = main_eblocks,
            .alias = "eblocks",
            .help = "Find maximal haplotype blocks."
        },
        {
            .func = main_ghmap,
            .alias = "ghmap",
            .help = "Haplotype association testing (ANOVA)."
        },
        {
            .func = main_convert,
            .alias = "convert",
            .help = "Convert VCF to NIEHS compact format with filtering (eblocks input file)"
        },
        {
            .func = main_annot,
            .alias = "annotate",
            .help = "Generate annotation file for variants from ensemble-VEP output file."
        },
        {
            .func = main_eigen,
            .alias = "pca",
            .help = "Principal component analysis."
        },
        {
            .func=NULL,
            .alias=NULL,
            .help=NULL
        }
};

static void usage()
{
    int i = 0;
    const char * sep = NULL;
    while (cmds[i].alias)
    {
        if (!cmds[i].func) sep = cmds[i].alias;
        if (sep)
        {
            std::cout<<sep<<" : "<<cmds[i].help<<std::endl;
        }
        i++;
    }
    std::cout << "Program: haplomap (haplotype-based computational genetic mapping, a.k.a HBCGM)\n"
    <<"Version: "<<__HAPLOMAP_VER__<<"\n\n"<<
    "Usage:   haplomap <subcommand> [options]\n\n"
    "Subcommands:\n\n"
    "    convert        convert and filter vcf (for eblocks -a)\n"
    "    annotate       annotate variants with ensembl-VEP result for (eblocks -g)\n"
    "    eblocks        find all maximal haploblocks\n"
    "    ghmap          haplotype association test (ANOVA)\n"
    "    pca            principal component analysis\n\n"
    "\nOptional arguments:\n"
    "    -v, --version  show program's version number and exit\n"
    "    -h, --help     show help message and exit."<< std::endl;
}

int main(int argc, char **argv) {

    if (argc < 2) { usage(); return 1; }

    if (std::strcmp(argv[1], "version") == 0 || std::strcmp(argv[1], "--version") == 0 || std::strcmp(argv[1], "-v") == 0)
    {
        std::cout <<
        "Program: haplomap (haplotype-based computational genetic mapping, a.k.a HBCGM)\n"
        "Version: "<<__HAPLOMAP_VER__<<"\n\n"<<
        "Compiled by "<<__COMPILER__<<" "<<__VERSION__<<std::endl;
        return 0;
    }
    else if (std::strcmp(argv[1], "help") == 0 || std::strcmp(argv[1], "--help") == 0 || std::strcmp(argv[1], "-h") == 0)
    {
        if (argc == 2) { usage(); return 0; }
        // Otherwise change "haplomap help COMMAND [...]" to "haplomap COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }

    int i = 0;
    while (cmds[i].alias)
    {
        if (cmds[i].func && strcmp(argv[1], cmds[i].alias)==0)
        {
            return cmds[i].func(argc-1,argv+1);
        }
        i++;
    }
    std::cerr<<"Error:: "<<__func__<<" unrecognized command: "<<argv[1]<<std::endl;
    return 1;

}
