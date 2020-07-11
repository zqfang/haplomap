

#include <getopt.h>
#include <cstring>
#include <unordered_map>
#include "ghmap.h"
#include "manova.h"
#include "fdr.h"
#include "haplomap.h"

struct GhmapOptions
{
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
    char *geneticRelationMatrix;
    char *geneticRelationIDs;
    // constructor
    GhmapOptions() : isCategorical(false), filterCoding(false), haploBlocks(false),
                     geneBlocks(false), pvalueCutoff(0.05),datasetName((char *)"Unnamed_dataset"),
                     phenotypeFileName(NULL), blocksFileName(NULL), outputFileName(NULL), geneName(NULL),
                     equalFile(NULL), goTermFile(NULL), goFilter(NULL),
                     geneticRelationMatrix(NULL), geneticRelationIDs(NULL) {};
};

GhmapOptions *parseGhmapOptions(int argc, char **argv)
{
    int c;

    GhmapOptions *opts = new GhmapOptions();

    static struct option long_options_ghmap[] = {
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
            {"relation", required_argument, 0, 'r'},
            /*{"version", no_argument, 0, 0},*/
            {0, 0, 0, 0}};

    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "cfkmahvn:p:b:l:o:g:e:q:t:i:r:", long_options_ghmap, &option_index);

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
                std::cout << "usage: ghmap [options]\n"
                        "\nrequired arguments:\n"  
                        "    -p, --phenotypes_file  <phenotype file> strain order should match to (-b)\n"
                        "    -b, --blocks_file      <output file from eblocks>\n"
                        "    -r, --relation         <genetic relation file>  n x n matrix\n"
                        "    -o, --output_file      <output file name>\n"
                        "\noptional arguments:\n"     
                        "    -n, --name             <name of phenotype dataset>\n"  
                        "    -e, --expression_file  <name of file>\n" 
                        "    -q, --equal_file       <name of file>\n"
                        "    -t, --goterms_file     <name of file>\n"
                        "    -g, --gene             <name of file>  writing block-oriented results file for gene\n"
                        "    -i, --goterms_include_file <name of file>\n"   
                        "                           output only genes with these terms\n"
                        "    -c    phenotype (-p) is categorical\n"                        
                        "    -f    filter out non-coding blocks\n"
                        "    -k    haploblocks generate a blocks-oriented results\n"
                        "    -m    output gene and haplotype block\n"
                        "    -a    output gene and haplotype block sort by block\n"
                        "    -l, --pvalue_cutoff    only write results with pvalue < cutoff\n"
                        "    -v, --verbose\n" 
                        "    -h, --help\n"
                     << std::endl;
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
            case 'r':
            {
                opts->geneticRelationMatrix = optarg;
                ///MARK: if set .id file here, cmdline could not parse relative path
                //char* ids = strdup(optarg);
                //opts->geneticRelationIDs = std::strcat(ids, ".id");
                break;
            }
            case '?':
            {
                /* getopt_long already printed an error message. */
                break;
            }
//            case 0: /* long option without a short arg */
//            {
//                if (strcmp("version", long_options_ghmap[option_index].name) == 0)
//                {
//                    std::cout << __COMPILER__ <<" "<< __VERSION__ << std::endl;
//                    std::cout << "ghmap version: " << __HAPLOMAP_VER__ << std::endl;
//                }
//                exit(1);
//            }
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
    } else if (NULL == opts->geneticRelationMatrix)// || opts->geneticRelationIDs == NULL)
    {
        std::cout<<"Required arg missing: genetic relation files (.rel, .rel.id)" << std::endl;
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

int main_ghmap(int argc, char **argv)
{

    GhmapOptions *opts = parseGhmapOptions(argc, argv);

    beginPhase("reading blocks summary file");
    readBlockSummary(opts->blocksFileName, opts->geneName, opts->goTermFile);
    endPhase();

    std::vector<std::vector<float>> phenvec(numStrains);
    beginPhase("reading genetic relation matrix");
    char* _relid = strdup(opts->geneticRelationMatrix);
    opts->geneticRelationIDs = std::strcat(_relid, ".id");
    //const_cast<char*>(mid.c_str()); // string to char*
    //printf("rel: %s\n",opts->geneticRelationMatrix);
    //printf("rel.id: %s\n",opts->geneticRelationIDs);

    MANOVA aov(opts->geneticRelationMatrix, opts->geneticRelationIDs, 4);
    aov.setEigen(); // calculate eigenvectors
    endPhase();

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
      bool ok = aov.setNonQMarkMat(pBlock->pattern, strainAbbrevs);
      if (ok)
          aov.pillaiTrace(pBlock->relFStat, pBlock->relPvalue);

      if (pBlock->FStat == INFINITY && pBlock->effect < 0.0)
      {
        cout << "Weird effect:" << endl;
        cout << "blockIdx = " << pBlock->blockIdx << ", FStat = " << pBlock->FStat << ", effect = " << pBlock->effect << endl;
      }
    }
    }
    endPhase();
    setBlockStats();
    beginPhase("Benjamini Hochberg procedure for controlling the FDR");
    bh_fdr(blocks, 0.05, 0);
    if (!opts->isCategorical)
        bh_fdr(blocks, 0.05,1);
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
