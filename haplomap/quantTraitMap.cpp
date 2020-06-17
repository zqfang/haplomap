

#include <getopt.h>
#include <cstring>
#include <unordered_map>
#include "ghmap.h"

#ifndef __HAPLOMAP_VER__
#define __HAPLOMAP_VER__ "0.1.0"
#endif

#ifdef __clang__
#define __COMPILER__ "clang++"
#else
#define __COMPILER__ "g++"
#endif


class Options
{

public:
  bool isCategorical;
  bool filterCoding;
  bool haploBlocks;
  bool geneBlocks;
  bool geneByBlocks;
  float pvalueCutoff;
  char *datasetName;
  char *phenotypeFileName;
  char *blocksFileName;
  char *outputFileName;
  char *geneName;
  char *expressionFile;
  char *equalFile;
  char *goTermFile;
  char *goFilter;
  // constructor
  Options() : isCategorical(false), filterCoding(false), haploBlocks(false), geneBlocks(false), pvalueCutoff(0.05),
              datasetName((char *)"Unnamed_dataset"), phenotypeFileName(NULL),
              blocksFileName(NULL), outputFileName(NULL), geneName(NULL),
              equalFile(NULL), goTermFile(NULL), goFilter(NULL){};
};

// Options parsing
Options *parseOptions(int argc, char **argv)
{
  int c;

  Options *opts = new Options();

  static struct option long_options[] = {
      // I wanted to use the flag feature, but it wouldn't work.
      {"help", no_argument, 0, 'h'},
      {"verbose", no_argument, 0, 'v'},
      {"categorical", no_argument, 0, 'c'},
      {"filter_coding", no_argument, 0, 'f'},
      {"haploblocks", no_argument, 0, 'k'},
      {"gene_block", no_argument, 0, 'm'},
      {"gene_block_by_block", no_argument, 0, 'a'},
      {"pvalue_cutoff", no_argument, 0, 'l'},
      {"name", required_argument, 0, 'n'},
      {"phenotypes_file", required_argument, 0, 'p'},
      {"blocks_file", required_argument, 0, 'b'},
      {"output_file", required_argument, 0, 'o'},
      {"gene", required_argument, 0, 'g'},
      {"expression_file", required_argument, 0, 'e'},
      {"equal_file", required_argument, 0, 'q'},
      {"goterms_file", required_argument, 0, 'q'},
      {"goterms_include_file", required_argument, 0, 'q'},
      {"version", no_argument, 0, 0},
      {0, 0, 0, 0}};

  while (1)
  {

    int option_index = 0;
    c = getopt_long(argc, argv, "cfkmahvn:p:b:l:o:g:e:q:t:i:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
    {
      break;
    }

    switch (c)
    {

    case 'c':
    {
      opts->isCategorical = true;
      break;
    }

    case 'f':
    {
      opts->filterCoding = true;
      break;
    }

    case 'k':
    {
      opts->haploBlocks = true;
      break;
    }

    case 'm':
    {
      opts->geneBlocks = true;
      break;
    }

    case 'a':
    {
      opts->geneByBlocks = true;
      break;
    }

    case 'l':
    {
      opts->pvalueCutoff = atof(optarg);
      break;
    }

    case 'h':
    {
      cout << "usage:\n ghmap [-c]\n"
              "    -v --verbose\n"
              "    -n --name <name of dataset>\n"
              "    -f --filter out non-coding blocks\n"
              "    -k --haploblocks generate a blocks-oriented results\n"
              "    -m --output gene and haplotype block\n"
              "    -a --output gene and haplotype block sort by block\n"
              "    -p --phenotypes_file <file with phenotype data>\n"
              "    -o --output_file <output file name>\n"
              "    -c --categorical\n"
              "    -l --pvalue_cutoff\n"
              "    -b --blocks_file\n"
              "    -g --gene writing block-oriented results file for gene\n"
              "    -e --expression_file\n"
              "    -q --equal_file <name of file>\n"
              "    -t --goterms_file <name of file>\n"
              "    -i --goterms_include_file output only genes with these terms <name of file>\n"
              "    --version\n"
           << endl;
      break;
    }

    case 'v':
    {
      verbose = true;
      break;
    }

    case 'n':
    {
      // cout << "option -n with arg " << optarg << endl;
      opts->datasetName = optarg;
      break;
    }

    case 'o':
    {
      // cout << "option -o with arg " << optarg << endl;
      opts->outputFileName = optarg;
      break;
    }

    case 'p':
    {
      // cout << "option -p with arg " << optarg << endl;
      opts->phenotypeFileName = optarg;
      break;
    }

    case 'b':
    {
      // cout << "option -b with arg " << optarg << endl;
      opts->blocksFileName = optarg;
      break;
    }

    case 'g':
    {
      // cout << "option -g with arg " << optarg << endl;
      opts->geneName = optarg;
      break;
    }

    case 'e':
    {
      // cout << "option -e with arg " << optarg << endl;
      opts->expressionFile = optarg;
      break;
    }

    case 'q':
    {
      // cout << "option -q with arg " << optarg << endl;
      opts->equalFile = optarg;
      break;
    }
    case 't':
    {
      opts->goTermFile = optarg;
      break;
    }
    case 'i':
    {
      opts->goFilter = optarg;
      break;
    }
    case '?':
    {
      /* getopt_long already printed an error message. */
      break;
    }
    case 0: /* long option without a short arg */
    {
      if (strcmp("version", long_options[option_index].name) == 0)
      {
        std::cout << __COMPILER__ <<" "<< __VERSION__ << std::endl;
        std::cout << "ghmap version: " << __HAPLOMAP_VER__ << std::endl;
      }
      exit(1);
    }
    default:
      abort();
    }
  }

  if (NULL == opts->blocksFileName)
  {
    cout << "Required arg missing: blocks file name (-b)" << endl;
    exit(1);
  }
  else if (NULL == opts->phenotypeFileName)
  {
    cout << "Required arg missing: phenotype file name (-p)" << endl;
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

int main(int argc, char **argv)
{
    // defined constants
//    const char *HAPLOCOLORS[] = {"red", "blue", "green", "orange", "violet", "yellow"};
//    int numCategories = 1; // default for non-categorical data.
//    vector<string> catNames; // maps strIdx -> category name.
//    // Maps gene names to a string of A's, M's, and P's
//    unordered_map<string, string> geneExprMap(40000);
//    // Globals
//    std::unordered_map<std::string, GeneSummary *> geneTable; // for gene-oriented interface
//    std::vector<BlockSummary *> blocks; // global vector of all blocks.

  Options *opts = parseOptions(argc, argv);

  beginPhase("reading blocks summary file");
  readBlockSummary(opts->blocksFileName, opts->geneName, opts->goTermFile);
  endPhase();


//BlockSummary *pBtest = blocks[1];
//cout << "Block: " << pBtest->chrName << "\t" << pBtest->blockIdx << endl;


  vector<vector<float>> phenvec(numStrains);

  beginPhase("reading phenotype file");
  if (opts->isCategorical)
  {
    readCPhenotypes(opts->phenotypeFileName, phenvec);
  }
  else
  {
    readQPhenotypes(opts->phenotypeFileName, phenvec);
  }
  endPhase();

  beginPhase("reading gene expressions file");
  if (opts->expressionFile)
  {
    readCompactGeneExpr(opts->expressionFile);
  }
  endPhase();

  beginPhase("reading go term file");
  if (opts->goTermFile)
  {
    vector<string> goTerms;
    readFileToVec(opts->goFilter, goTerms);
    filterGoTerms(opts->goTermFile, goTerms);
  }
  endPhase();

  if (opts->filterCoding)
  {
    beginPhase("filtering blocks by coding");
    filterCodingBlocks();
    endPhase();
  }

  if (opts->equalFile)
  {
    beginPhase("filtering by equality class");
    vector<int> equalClass;
    readEqualFile(opts->equalFile, equalClass);
    filterEqualBlocks(equalClass);
    endPhase();
  }

  beginPhase("computing ANOVA p-values");
  for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
  {
    BlockSummary *pBlock = blocks[blkIdx];
    //cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << endl;
    if (pBlock->isIgnored)
    {
      if (traceFStat)
      {
        cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << "\tINGORED" << endl;
      }
    }
    else
    {
      if (traceFStat)
      {
        cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << "\t";
        showPattern(pBlock->pattern);
        cout << endl;
      }
      ANOVA(phenvec, pBlock->pattern, pBlock->FStat, pBlock->pvalue, pBlock->effect);
      if (pBlock->FStat == INFINITY && pBlock->effect < 0.0)
      {
        cout << "Weird effect:" << endl;
        cout << "blockIdx = " << pBlock->blockIdx << ", FStat = " << pBlock->FStat << ", effect = " << pBlock->effect << endl;
      }
    }
  }
  endPhase();
  setBlockStats();
  beginPhase("setting block stats");

  endPhase();

  beginPhase("sorting blocks");
  BlocksComparator bcomp(opts->isCategorical);
  sort(blocks.begin(), blocks.end(), bcomp);
  endPhase();

  beginPhase("sorting blocks in gene table");
  // Pass to sort block vectors in the gene table.
  for (std::unordered_map<string, GeneSummary *>::iterator git = geneTable.begin(); git != geneTable.end(); git++)
  {
    vector<BlockSummary *> &gBlocks = (*git).second->blocks;
    sort(gBlocks.begin(), gBlocks.end(), bcomp);
  }
  endPhase();

  // If haploBlocks flag is set, generate a blocks-oriented results file instead a gene-oriented file.
  // This is a bit of a hack.  If geneName is provided on command line, this generates a block-oriented
  // results file, in the format of nhaplomap.pl, for display when someone clicks on a gene in the
  // gene-oriented html (or if * was specified for gene name)
  // Otherwise, it generates the gene-oriented html.
  if (opts->geneName || opts->haploBlocks)
  {
    beginPhase("writing block-oriented results file for gene.");
    if (opts->isCategorical)
    {
      // Find FStat cutoff
      int cutoffBlockIdx = (int)(opts->pvalueCutoff * blocks.size());
      float FCutoff = blocks[cutoffBlockIdx]->FStat;
      //    cout << "pvalueCutoff = " << opts->pvalueCutoff << ", cutoffBlockIdx = " << cutoffBlockIdx
      //	 << ", num blocks = " << blocks.size()
      //	 << ", FCutoff = " << FCutoff << endl;
      writeBlockSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, FCutoff);
    }
    else
    {
      writeBlockSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, opts->pvalueCutoff);
    }
  }
  else if (opts->geneBlocks)
  {
    writeGeneBlockSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, opts->pvalueCutoff);
  }
  else if (opts->geneByBlocks)
  {
    writeGeneBlockByBlocks(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, opts->pvalueCutoff);
  }
  else
  {
    beginPhase("writing gene-oriented results file.");
    if (opts->isCategorical)
    {
      // Find FStat cutoff
      int cutoffBlockIdx = (int)(opts->pvalueCutoff * blocks.size());
      float FCutoff = blocks[cutoffBlockIdx]->FStat;
      writeGeneSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, FCutoff, opts->filterCoding);
    }
    else
    {
      writeGeneSums(opts->isCategorical, opts->outputFileName, opts->datasetName,
                    phenvec, blocks, opts->pvalueCutoff, opts->filterCoding);
    }
  }

  endPhase();

  return 0;
}
