//
// Created by Zhuoqing Fang on 7/11/20.
//

#ifndef __HAPLOMAP_VER__
#define __HAPLOMAP_VER__ "0.1.0"
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
            .help = "find haploblocks"
        },
        {
            .func = main_ghmap,
            .alias = "ghmap",
            .help = "run ghmap"
        },
        {
            .func=NULL,
            .alias=NULL,
            .help=NULL
        }
};

static void usage(){

    int i = 0;
    const char * sep = NULL;
    while (cmds[i].alias)
    {
        if (!cmds[i].func) sep = cmds[i].alias;
        if (sep)
        {
            std::cout<<"\n"<<sep<<std::endl;
        }
        i++;
    }
    std::cout <<
    "Program: haplomap (haplotype-based computational genetic mapping, a.k.a HBCGM)\n"
    "Version: "<<__HAPLOMAP_VER__<<"\n"<<
    "Compiled by "<<__COMPILER__<<" "<<__VERSION__<<std::endl;
    std::cout <<"\n"<<
    "Usage:    haplomap <subcommand> [options]\n\n"
    "Subcommands:\n"
    "    eblocks        find all maximal haploblocks\n"
    "    ghmap          statistical testing for hapbloblock with trait\n"
    "Optional arguments:\n"
    "    -v, --version  show program's version number and exit\n"
    "    -h --help      show help message and exit."
              << std::endl;
}

int main(int argc, char **argv) {

    if (argc < 2) { usage(); return 1; }

    if (std::strcmp(argv[1], "version") == 0 || std::strcmp(argv[1], "--version") == 0 || std::strcmp(argv[1], "-v") == 0)
    {
        std::cout<< "HAPLOMAP: "<<  __HAPLOMAP_VER__ << std::endl;
        std::cout<<__COMPILER__<<" "<< __VERSION__<<std::endl;
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
