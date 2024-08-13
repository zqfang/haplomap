

#include <getopt.h>
#include <cstring>
#include <unordered_map>
#include <memory>
#include "ghmap.h"
#include "stats.h"
#include "haplomap.h"

struct GhmapOptions
{
    bool isCategorical;
    bool filterCoding;
    bool haploBlocks;
    bool geneBlocks;
    bool geneAllBlocks;
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
    // constructor
    GhmapOptions() : isCategorical(false), filterCoding(false), haploBlocks(false),
                     geneBlocks(false), geneAllBlocks(false), pvalueCutoff(0.05),datasetName((char *)"Unnamed_dataset"),
                     phenotypeFileName(NULL), blocksFileName(NULL), outputFileName(NULL), geneName(NULL),
                     expressionFile(NULL), equalFile(NULL), goTermFile(NULL), goFilter(NULL),
                     geneticRelationMatrix(NULL) {};
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
            {"gene_all_blocks", no_argument, 0, 'a'},
            {"pvalue_cutoff", no_argument, 0, 'l'},// optional_argument means without a value is OK.
            {"name", required_argument, 0, 'n'}, // required_argument means an associated value must given
            {"phenotypes", required_argument, 0, 'p'},
            {"blocks", required_argument, 0, 'b'},
            {"output", required_argument, 0, 'o'},
            {"gene", required_argument, 0, 'g'},
            {"expression", required_argument, 0, 'e'},
            {"equal", required_argument, 0, 'q'},
            {"goterms", required_argument, 0, 't'},
            {"goterms_include", required_argument, 0, 'i'},
            {"relation", required_argument, 0, 'r'},
            /*{"version", no_argument, 0, 0},*/
            {0, 0, 0, 0}};

    const char *usage = "usage: ghmap [options]\n"
                        "\nRequired arguments:\n"
                        "    -p, --phenotypes       <phenotype file> \n"
                        "                                The same input file of (eblocks -b)\n"
                        "    -b, --blocks           <haploblocks file >\n"
                        "                                The output file from (eblocks -o)\n"
                        "    -o, --output           <output file name>\n"
                        "                                Output gene-summaried results by default.\n"
                        "\nOptional arguments:\n"
                        "    -n, --name             <name of phenotype dataset> \n"
                        "                                 To include INDEL/SV codon flag, please make sure to \n" 
                        "                                 add _INDEL or _SV as a suffix to the name explicitly. \n"
                        "                                 e.g. MPD123_SV. Default, output SNP CodonFlag .\n"
                        "    -c, --categorical           phenotype (-p) is categorical\n"
                        "    -r, --relation         <genetic relation file> \n" 
                        "                                <.rel> file for population structure analysis.\n"
                        "                                n x n matrix with header line (startswith '#') contain sample names.\n"
                        "    -e, --expression       <name of file>\n"
                        "    -q, --equal            <name of file>\n"
                        "    -t, --goterms          <name of file>\n"
                        "    -i, --goterms_include  <name of file> \n"
                        "                                 Output only genes with these terms\n"
                        "    -f, --filter_coding          Filter out non-coding blocks\n"
                        "    -g, --gene                   Output gene-summaried results. Default.\n"
                        "                                 NOTE:: Only write the overlapped halpoblock with best pvalue/Fstat, representing all overlapped blocks. \n"
                        "                                     The CodonFlag is an aggregated indicator showing that a gene has blocks with coding change. \n"
                        "                                     The best block itself might not contain any coding changes. \n"
                        "                                     Run the ghmap with -a/k/m tag will give you all overlapped blocks with correct CodonFlag. \n"
                        "    -a, --gene_all_blocks       Output gene-oriented results of all blocks that overalp a gene. \n"
                        "    -m, --gene_block            Output gene-oriented results block by block. almost the same to -a\n"
                        "    -k, --haploblocks           Output block-oriented results.\n"
                        "    -l, --pvalue_cutoff    <float>   Only write results with pvalue < cutoff. Default: 0.05\n"
                        "    -v, --verbose\n"
                        "    -h, --help\n";

    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "hvcfkman:p:b:l:o:g:e:q:t:i:r:", long_options_ghmap, &option_index);

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
                opts->geneAllBlocks = true;
                break;
            }

            case 'l':
            {
                opts->pvalueCutoff = atof(optarg);
                break;
            }

            case 'h':
            {
                std::cout << usage << std::endl;
                std::exit(0);
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
                break;
            }
            case '?':
            {
                /* getopt_long already printed an error message. */
                break;
            }
            default:
                std::abort();
        }
    }
    if (argc == 1)
    {
        std::cout<<usage<<std::endl;
        std::exit(1);
    }

    if (NULL == opts->blocksFileName)
    {
        std::cout<<usage<<std::endl;
        std::cout << "Required arg missing: blocks file name (-b)" << std::endl;
        std::exit(1);
    }
    else if (NULL == opts->phenotypeFileName)
    {
        std::cout<<usage<<std::endl;
        std::cout << "Required arg missing: phenotype file name (-p)" << std::endl;
        std::exit(1);
    }
//    else if (NULL == opts->geneticRelationMatrix)
//    {
//        std::cout<<usage<<std::endl;
//        std::cout<<"Required arg missing: genetic relation files (.rel, .rel.id) -r "<< std::endl;
//        exit(1);
//    }

    // Print any remaining command line arguments (not options).
    if (optind < argc)
    {
        std::cout << "Extraneous things on command line: ";
        while (optind < argc)
        {
            std::cout << argv[optind++] << std::endl;
        }
        exit(1);
    }

    return opts;
}

int main_ghmap(int argc, char **argv)
{

    GhmapOptions *opts = parseGhmapOptions(argc, argv);

    beginPhase("reading blocks summary file");
    readBlockSummary(opts->blocksFileName); // store in blocks
    endPhase();
    /// global variable defined in haplolib
    // Dynum<std::string> strainAbbrevs;
    // int numStrains = -1;
    std::vector<std::vector<float>> phenvec(numStrains);
    beginPhase("reading phenotype file");
    if (opts->isCategorical)
    {
        readCPhenotypes(opts->phenotypeFileName, phenvec, strainAbbrevs);
    }
    else
    {
        readQPhenotypes(opts->phenotypeFileName, phenvec, strainAbbrevs);
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
        std::vector<std::string> goTerms;
        readFileToVec(opts->goFilter, goTerms);
        filterGoTerms(opts->goTermFile, goTerms);
    }
    endPhase();


    if (opts->equalFile)
    {
        beginPhase("filtering by equality class");
        std::vector<int> equalClass;
        readEqualFile(opts->equalFile, equalClass);
        filterEqualBlocks(equalClass);
        endPhase();
    }

    std::shared_ptr<MANOVA> manova;
    std::shared_ptr<ANOVA> anova;
    anova = std::make_shared<ANOVA>(phenvec);

    if (opts->geneticRelationMatrix != NULL) {
        beginPhase("reading genetic relation matrix");
        manova = std::make_shared<MANOVA>(opts->geneticRelationMatrix, (char *)"\t",  4);
        manova->setEigen(strainAbbrevs); // calculate eigenvectors for selected strains
        endPhase();
    }

    beginPhase("computing ANOVA p-values");
    for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
    {
        BlockSummary *pBlock = blocks[blkIdx];
        // ANOVA analysis
        anova->stat(pBlock->pattern, pBlock->FStat, pBlock->pvalue, pBlock->effect);
        // population structure analysis
        if (opts->geneticRelationMatrix != NULL) 
        {
            bool ok = manova->setNonQMarkMat(pBlock->pattern, strainAbbrevs);
            if (ok)
                manova->pillaiTrace(pBlock->relFStat, pBlock->relPvalue);
            pBlock->relIgnore = false; // means final output will contain pop analysis results
        }
        if (pBlock->FStat == INFINITY && pBlock->effect < 0.0)
        {
            cout << "Weird effect:" << endl;
            cout << "blockIdx = " << pBlock->blockIdx << ", FStat = " << pBlock->FStat << ", effect = " << pBlock->effect << endl;
        }
        //cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << endl;
        if (pBlock->isIgnored && traceFStat)
        {
            std::cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << "\t";
            pBlock->showPatten();
            std::cout << endl;
        
        }
    }
    endPhase();

    beginPhase("Benjamini Hochberg procedure for controlling the FDR");
    if (opts->geneticRelationMatrix != NULL)
        bh_fdr(blocks, 0.05, 0); // pop structure pval correction
    if (!opts->isCategorical)
        bh_fdr(blocks, 0.05, 1); // annova pval correction
    endPhase();

    beginPhase("sorting blocks");
    BlocksComparator bcomp(opts->isCategorical);
    std::sort(blocks.begin(), blocks.end(), bcomp);
    endPhase();

    beginPhase("sorting blocks in gene table");
    // Pass to sort block vectors in the gene table.
    for (std::unordered_map<string, GeneSummary *>::iterator git = geneTable.begin(); git != geneTable.end(); git++)
    {
        vector<BlockSummary *> &gBlocks = (*git).second->blocks;
        std::sort(gBlocks.begin(), gBlocks.end(), bcomp);
    }
    endPhase();

    // If haploBlocks flag is set, generate a blocks-oriented results file instead a gene-oriented file.
    // This is a bit of a hack.  If geneName is provided on command line, this generates a block-oriented
    // results file, in the format of nhaplomap.pl, for display when someone clicks on a gene in the
    // gene-oriented html (or if * was specified for gene name)
    // Otherwise, it generates the gene-oriented html.
    float cutoff = opts->pvalueCutoff;
    if (opts->isCategorical)
    {
        // Find FStat cutoff
        int cutoffBlockIdx = (int)(opts->pvalueCutoff * blocks.size());
        cutoff = blocks[cutoffBlockIdx]->FStat;
    }

    // if ( opts->geneName || opts->haploBlocks) // -k
    if ( opts->haploBlocks) // -k , click on a gene in nhaplomap.pl will need to trigger -k now 
    {
        beginPhase("writing block-oriented results file for gene.");
        writeBlockSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, cutoff);
    }
    else if (opts->geneBlocks) // -m
    {
        beginPhase("writing gene-oriented results file for gene.");
        writeGeneBlockSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, cutoff);
    }
    else if (opts->geneAllBlocks)  // -a
    {
        beginPhase("writing gene-oriented results file.");
        writeGeneBlockByBlocks(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, cutoff, opts->filterCoding);
    }
    else // if (opts->geneName)
    {
        beginPhase("writing gene-summaried results file.");
        writeGeneSums(opts->isCategorical, opts->outputFileName, opts->datasetName, phenvec, blocks, cutoff, opts->filterCoding);
    }

    endPhase();

    // free memory
    delete opts;
    for (auto & block : blocks)
    {
        delete block;
    }

    for (auto & git : geneTable){
        delete git.second;
    }

    return 0;
}
