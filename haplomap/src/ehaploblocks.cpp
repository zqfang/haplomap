// -*-c++-*-

// (c) 2008, 2009, 2010 Board of Trustees of Leland Stanford Jr. University.  All rights reserved.
// Author: David L. Dill

// Fast algorithm to find haplotype blocks ala Peltz.

// ISSUES:
//    SNPs are in somewhat random order w/in chromosome
//    Gene names are repeated for a single snp (Oprm1)
//    Gene names are in various orders.

// TODO:
// FIXME: diplay blocks in a specified segment
// FIXME: non-monotonicity bug

// This version of the code changed the block finding significantly.  It finds all maximum-length blocks
// for 2,3,4,5 haplotypes from a given position, and doesn't worry about scores.
#include <getopt.h>
#include <unordered_map>
#include <cstring>
#include "eblocks.h"
#include "haplomap.h"


struct EblockOptions
{
    bool non_overlapping;
    char *strainsFileName;
    char *allelesFileName;
    char *genesFileName;
    char *blocksFileName;
    char *blockSNPsFileName;
    int minBlockSNPs;

    // constructor
    EblockOptions() : strainsFileName(NULL), blocksFileName(NULL){};
};



EblockOptions *parseEblockOptions(int argc, char **argv)
{
    int c;

    EblockOptions *opts = new EblockOptions();

    static struct option long_options_eblock[] = {
            // I wanted to use the flag feature, but it wouldn't work.
            {"help", no_argument, 0, 'h'},
            {"verbose", no_argument, 0, 'v'},
            {"non_overlapping", no_argument, 0, 'n'},
            {"strains_file", required_argument, 0, 's'},
            {"alleles_file", required_argument, 0, 'a'},
            {"genes_file", required_argument, 0, 'g'},
            {"blocks_file", required_argument, 0, 'o'},
            {"blockSNPs_file", required_argument, 0, 'p'},
            {"min_block_SNPs", required_argument, 0, 'm'},
            /*{"version", no_argument, 0, 0},*/
            {0, 0, 0, 0}};

    // default values
    opts->minBlockSNPs = 4;

    const char *usage = "usage: eblocks [options]\n"
                        "\nrequired arguments:\n" 
                        "    -s, --strains_file     <name of strains file>\n"
                        "    -a, --alleles_file     <name of alleles file>\n"
                        "    -g, --genes_file       <name of SNPs-to-genes file>\n"
                        "    -o, --blocks_file      <output file with haplotype blocks>\n"

                        "\noptional arguments:\n"                       
                        "    -p, --blockSNPs_file   <output file with SNPs for blocks>\n"
                        "    -m, --min_block_SNPs\n"
                        "    -n, --non_overlapping\n" 
                        "    -v, --verbose\n"
                        "    -h, --help\n";

    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "hvna:g:s:o:p:m:", long_options_eblock, &option_index);

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
                verbose = true;
                break;
            }

            case 'n':
            {
                opts->non_overlapping = true;
                break;
            }

            case 's':
            {
                opts->strainsFileName = optarg;
                break;
            }

            case 'a':
            {
                opts->allelesFileName = optarg;
                break;
            }

            case 'g':
            {
                opts->genesFileName = optarg;
                break;
            }

            case 'o':
            {
                opts->blocksFileName = optarg;
                break;
            }

            case 'p':
            {
                opts->blockSNPsFileName = optarg;
                break;
            }

            case 'm':
            {
                // cout << "option -m with arg " << optarg << endl;
                opts->minBlockSNPs = atoi(optarg);
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

    // isn't there some way to do this automatically based on above info?
    if (NULL == opts->strainsFileName)
    {
        cerr << usage << endl;
        cerr << "Required arg missing: strains file name (-s)" << endl;
        exit(1);
    }
    else if (NULL == opts->allelesFileName)
    {
        cerr << usage << endl;
        cerr << "Required arg missing: alleles file name (-a)" << endl;
        exit(1);
    }
    else if (NULL == opts->blocksFileName)
    {
        cerr << usage << endl;
        cerr << "Required arg missing: blocks file name (-o)" << endl;
        exit(1);
    }
    // else if (NULL == opts->blockSNPsFileName)
    // {
    //     cerr << "Required arg missing: block SNPs file name (-p)" << endl;
    //     cerr << usage << endl;
    //     exit(1);
    // }

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


// Main for nhaploblocks
int main_eblocks(int argc, char **argv)
{
    EblockOptions *opts = parseEblockOptions(argc, argv);

  beginPhase("reading strains");
  readStrains(opts->strainsFileName);
  endPhase();
  numStrains = strainAbbrevs.size(); 
  // Can't do this earlier because numStrains is used in constructor.
  snpVec.reserve(approxNumSNPs);

  // FIXME: Add option for "compact SNP file"
  beginPhase("reading compact alleles file");
  readAlleleInfoCompact(opts->allelesFileName);
  endPhase();
  std::string chr = chromosomes.eltOf(0);
  // if no gene names file specified, don't read it, omit gene names.
  if (opts->genesFileName)
  {
    beginPhase("reading SNP gene names file");
    // readAlleleInfoCompact(opts->geneNamesFile);
    readSNPGeneNames(opts->genesFileName);
    endPhase();
  }

  minDefined = (numStrains + 1) / 2; // half of strains must be defined (rounded up).

  beginPhase("filtering and sorting SNPs");
  filterAndSortSNPs();
  endPhase();

  if (opts->non_overlapping)
  {
    beginPhase("finding compound blocks");
    findCompoundBlocks(haploLimit);
    endPhase();

    beginPhase("choosing non-overlapping blocks");
    chooseBlocks();
    endPhase();
  }
  else
  {
    // Current method (overlapping blocks ok).
    beginPhase("finding maximal blocks");
    findAllMaximalBlocks(haploLimit, opts->minBlockSNPs);
    endPhase();
  }

  // showBlocksSNPs(chosenHaploBlocks);
  beginPhase("writing block summary");
  writeBlockSummary(opts->blocksFileName, opts->minBlockSNPs);
  endPhase();

  if (opts->blockSNPsFileName != NULL)
  {
      beginPhase("writing block SNPs");
      writeBlockSNPs(opts->blockSNPsFileName, opts->minBlockSNPs);
      endPhase();
  }

  if (verbose)
  {
    updateStats();
    reportStats();
  }
  delete opts;
  for (auto & snp: snpMap)
      delete snp.second;
  for (auto & b: chosenHaploBlocks)
      delete b;

  return 0;
}

//// This just compacts the perlegen tables.
//// KEEP THIS -- needed to compact perlegen tables
//int main3(int argc, char **argv)
//{
//
//  //  cout << "Just doing chromosome 1, for testing" << endl;
//  //  beginPhase("reading chr_index.txt");
//  //  readDynumList("PERLEGEN/chr_index.txt", chromosomes);
//  //  endPhase();
//
//  beginPhase("reading strains");
//  // read all perlegen strains.
//  readStrains((char *)"all_strains.txt");
//  endPhase();
//
//  numStrains = strainAbbrevs.size();
//
//  // perlegen_b36_snp_vs_mmgene_091208.unl
//  beginPhase("reading perlegen_b36_snp_vs_mmgene_091208.unl");
//  readPerlegenSNPvsmgene((char *)"perlegen_b36_snp_vs_mmgene_091208.unl");
//  // readPerlegenSNPvsmgene("perlegen_test_chr.unl");
//  endPhase();
//
//  beginPhase("reading perlegen_b36_strain_091208.unl");
//  readPerlegenAlleleInfo((char *)"perlegen_b36_strain_091208.unl");
//  //readPerlegenAlleleInfo("perlegen_test_snps.unl");
//  endPhase();
//
//  beginPhase("writing perlegen_b36_compact.txt");
//  writeAlleleInfoCompact((char *)"perlegen_b36_compact.txt");
//  // readPerlegenAlleleInfo("perlegen_test_snps.unl");
//  endPhase();
//
//  return 0;
//}
//
//// compact the roche tables
//// KEEP THIS -- you need it to build the compact tables.
//int main4(int argc, char **argv)
//{
//
//  beginPhase("reading chr_index.txt");
//  // readDynumList("/home/dill/bin/data/chr_index.txt", chromosomes);
//  // readDynumList("/gourd/dill/bin/data/chr_index.txt", chromosomes);
//  readDynumList((char *)"chr_index.txt", chromosomes);
//  endPhase();
//  // this has the strains we actually want to use.
//  // readDynumList("strain_index.txt", relevantStrains);
//
//  //  beginPhase("reading h2_strain_index.txt");
//  // readStrains("h2_strain_index.txt");
//  beginPhase("reading strains file");
//  readStrains((char *)"roche_all_strains.txt");
//  endPhase();
//
//  numStrains = strainAbbrevs.size();
//  minDefined = (numStrains + 1) / 2; // half of strains must be defined (rounded up).
//
//  beginPhase("reading chr_info_perl.txt");
//  // readChromosomeInfo("/gourd/dill/bin/data/chr_info_perl.txt");
//  readChromosomeInfo((char *)"chr_info_perl.txt");
//  endPhase();
//
//  beginPhase("reading allele_info_perl.txt");
//  // readAlleleInfo("/gourd/dill/bin/data/allele_info_perl.txt");
//  readAlleleInfo((char *)"allele_info_perl.txt");
//  // readAlleleInfo("allele_info_test.txt");
//  endPhase();
//
//  beginPhase("writing roche_compact.txt");
//  writeAlleleInfoCompact((char *)"roche_compact.txt");
//  // readPerlegenAlleleInfo("perlegen_test_snps.unl");
//  endPhase();
//
//  return 0;
//}
