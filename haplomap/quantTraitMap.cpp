// -*- C++ -*-
#ifndef __HAPLOMAP_VER__
#define __HAPLOMAP_VER__ 0.1
#endif
// FIXME: Add tracing of block info.
#include <getopt.h>
#include <cstring>
#include <unordered_map>
#include "gsl/gsl_cdf.h"
#include "haploLib.h"

const char *HAPLOCOLORS[] = {"red", "blue", "green", "orange", "violet", "yellow"};

const int AACLASSES[] = {4, -1, 3, 1, 1, 4, 2, 0, 4, -1, 0, 4, 4, 2, -1, 4, 2, 0, 2, 2, -1, 4, 4, 5, 2, -1};

int numCategories = 1; // default for non-categorical data.

vector<string> catNames; // maps strIdx -> category name.

int numHaplotypes(char *pattern); // forward decl.
int interestingChanges(const map<string, string> &geneCodingMap);

// Maps gene names to a string of A's, M's, and P's
unordered_map<string, string> geneExprMap(40000);

void upcase(string &str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

// Read the file of gene expression data
void readCompactGeneExpr(char *fname)
{
  ColumnReader rdr(fname, (char *)"\t");
  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    // A typical line: "Myc	PAPAAAAAAAPAA"
    string geneName = rdr.getToken(0);
    string present = rdr.getToken(1);

    upcase(geneName);

    // Gene names in expr file are all upper case, so use upper case comparisons.
    // Don't really need to do upcase here, but doing it anyway in case someone
    // substitutes a file with mixed-case names.
    geneExprMap[geneName] = present;
  }
}

// renumber eqclasses in pattern so the increase from left to right
void writeSortedPattern(ostream &os, char *pattern, vector<int> &strOrderVec)
{
  int numHaplo = numHaplotypes(pattern);
  char *sortedEqMap = (char *)malloc(numHaplo);
  memset(sortedEqMap, '?', numHaplo);
  char curEq = 0;

  for (int strIdx = 0; strIdx < numStrains; strIdx++)
  {
    int str = strOrderVec[strIdx];
    int eqclass = (int)pattern[str];
    if ('?' != eqclass)
    {
      if ('?' == sortedEqMap[eqclass])
      {                                 // class hasn't been mapped yet.
        sortedEqMap[eqclass] = curEq++; // renumber for left-to-right color consistency
      }
      os << (char)(sortedEqMap[eqclass] + '0');
    }
    else
    {
      os << '?';
    }
  }
}

// summary of a block, read for the file.
class BlockSummary
{

public:
  char *chrName;
  int blockIdx;
  int blockStart;
  int blockSize;
  int chrBegin;
  int chrEnd;
  char *pattern;
  bool isIgnored; // for blocks that have been filtered out.
  float FStat;
  float pvalue;
  float effect;
  int numHaplo;
  int numInteresting;
  map<string, string> geneIsInteresting; // gene name -> codon change, by gene BY
  map<string, string> geneIsCodingMap;   // gene name -> coding bit
  // constructor
  BlockSummary(char *chrnm, int num, int start, int size,
               int chrbeg, int chrend, char *pat)
      : chrName(chrnm), blockIdx(num), blockStart(start), blockSize(size),
        chrBegin(chrbeg), chrEnd(chrend), pattern(pat), isIgnored(false),
        FStat(INFINITY), pvalue(1.0), effect(0.0), numHaplo(-1), numInteresting(-1){};

  string updateCodonScore(string str)
  {
    int count = 0;
    size_t pos = 0;
    size_t endpos = 0;
    while (true)
    {
      pos = str.find("<", pos + 1);
      if (pos == string::npos)
        break;
      endpos = str.find(">", pos + 1);
      int aa1 = str[pos - 1] - 'A';
      int aa2 = str[endpos + 5] - 'A';
      if (AACLASSES[aa1] != AACLASSES[aa2])
      {
        count++;
        int X = 'X' - 'A';
        if (aa1 == X || aa2 == X)
        { // stop codon
          return "stop_codon";
        }
      }
    }
    if (str.find("SPLICE_SITE") != string::npos)
    {
      return "splicing";
    }
    if (count > 0)
    {
      return "significant_codon_change";
    }
    return "non_significant_codon_change";
  }

  friend void updateGeneIsInteresting(BlockSummary *pb)
  { // BY
    for (map<string, string>::iterator giit = pb->geneIsCodingMap.begin(); giit != pb->geneIsCodingMap.end(); giit++)
    {
      if (giit->second == "0" || giit->second == "1")
      {
        pb->geneIsInteresting[giit->first] = "no_codon_change";
        continue;
      }
      pb->geneIsInteresting[giit->first] = pb->updateCodonScore(giit->second);
    }
  }

  // print a line of the blocks file.
  // blockIdx	blockStart	blockSize	chromosome	begin	end	pattern	pval	effect	genename genehascoding ...
  friend void showBlockSum(ostream &os, bool isCategorical, BlockSummary *pb, vector<int> &strOrderVec)
  {
    os << pb->blockIdx;
    if (pb->isIgnored)
    {
      os << "\t(IGNORED)";
    }
    os << "\t" << pb->blockStart << "\t" << pb->blockSize
       << "\t" << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t";

    writeSortedPattern(os, pb->pattern, strOrderVec);

    os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) << "\t" << pb->effect;
    updateGeneIsInteresting(pb);

    // gene names and coding bits
    map<string, string> &geneIsCodingMap = pb->geneIsCodingMap;
    for (map<string, string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
    { ///EDITED BY
      os << "\t" << (*giit).first << "\t" << (*giit).second << "\t" << pb->geneIsInteresting[(*giit).first];
    }
    os << endl;
  }
  // BY addition, output as such:
  // gene	codon_flag	pattern	pval	effect	chromosome	begin	end	blockIdx	blockStart	blockSize	expression
  friend void showGeneBlockByBlock(ostream &os, bool isCategorical, BlockSummary *pb, vector<int> &strOrderVec)
  {
    updateGeneIsInteresting(pb);
    map<string, string> &geneIsCodingMap = pb->geneIsCodingMap;
    map<string, string>::iterator giit = geneIsCodingMap.begin();
    if (giit == geneIsCodingMap.end())
    { // no genes in this block
      os << "None\tNone\t";
      writeSortedPattern(os, pb->pattern, strOrderVec);
      os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) << "\t" << pb->effect << "\t"
         << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t"
         << pb->blockIdx << pb->blockStart << "\t" << pb->blockSize << endl;
    }
    else
    {
      for (map<string, string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
      {
        os << (*giit).first << "\t" << pb->geneIsInteresting[(*giit).first] << "\t";
        writeSortedPattern(os, pb->pattern, strOrderVec);
        os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) << "\t" << pb->effect << "\t"
           << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t"
           << pb->blockIdx << pb->blockStart << "\t" << pb->blockSize;
        string gname = (*giit).first;
        upcase(gname);
        if (geneExprMap.find(gname) == geneExprMap.end())
        {
          os << "\t-------------" << endl;
        }
        else
        {
          os << "\t" << geneExprMap[gname] << endl;
        }
      }
    }
  }

  // Return true if pval is above cutoff or FStat is below it.
  friend bool isCutoff(bool isCategorical, float cutoff, BlockSummary *pBlock)
  {
    return (isCategorical) ? (pBlock->FStat < cutoff) : (pBlock->pvalue > cutoff);
  }

  friend void showBlockSums(ostream &os, bool isCategorical,
                            vector<BlockSummary *> &blocks, float cutoff, vector<int> &strOrderVec)
  {
    for (unsigned i = 0; i < blocks.size(); i++)
    {
      if (!isCutoff(isCategorical, cutoff, blocks[i]))
      {
        showBlockSum(os, isCategorical, blocks[i], strOrderVec);
      }
    }
  }

  friend void showGeneBlockByBlocks(ostream &os, bool isCategorical, vector<BlockSummary *> &blocks, float cutoff, vector<int> &strOrderVec)
  {
    for (unsigned i = 0; i < blocks.size(); i++)
    {
      if (!isCutoff(isCategorical, cutoff, blocks[i]))
      {
        showGeneBlockByBlock(os, isCategorical, blocks[i], strOrderVec);
      }
    }
  }
};

void showIsCoding(map<string, string> geneIsCodingMap)
{
  for (map<string, string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
  {
    cout << "\t" << (*giit).first << "\t" << (*giit).second;
  }
  cout << endl;
}

// This is to sort blocks in the best order.  First, sort by p-value
// or F statistic, depending on whether the phenotype is quantitative
// or categorical, then by chromosome, then by position.
class BlocksComparator
{
public:
  bool isCategorical;

  // constructor
  BlocksComparator(bool isCat) : isCategorical(isCat){};

  bool operator()(const BlockSummary *pb1, const BlockSummary *pb2) const
  {
    // Use F statistic or p value, depending on whether it's categorical or not.
    if (isCategorical)
    {
      if (pb1->FStat > pb2->FStat)
      {
        return true;
      }
      else if (pb1->FStat < pb2->FStat)
      {
        return false;
      }
    }
    else
    {
      if (pb1->pvalue < pb2->pvalue)
      {
        return true;
      }
      else if (pb1->pvalue > pb2->pvalue)
      {
        return false;
      }
    }

    // Break ties by choosing blocks with the fewest haplotypes.
    int nh1 = pb1->numHaplo; //numHaplotypes(pb1->pattern);
    int nh2 = pb2->numHaplo; //numHaplotypes(pb2->pattern);

    if (nh1 < nh2)
    {
      return true;
    }
    else if (nh2 < nh1)
    {
      return false;
    }

    // Then order by SNP changes
    int snpChangeScore1 = pb1->numInteresting; //interestingChanges(pb1->geneIsCodingMap);
    int snpChangeScore2 = pb2->numInteresting; //interestingChanges(pb2->geneIsCodingMap);

    if (snpChangeScore1 < snpChangeScore2)
    {
      return false;
    }
    else if (snpChangeScore2 < snpChangeScore1)
    {
      return true;
    }
    // Then order by chromosome
    char *chr1 = pb1->chrName;
    char *chr2 = pb2->chrName;

    if (strcmp(chr1, chr2) != 0)
    {
      // deal with complexities of comparing chr names.
      // put "M" last
      if (strcmp(chr2, "M") == 0)
      {
        return true;
      }
      else if (strcmp(chr1, "M") == 0)
      {
        return false;
      }
      // Y next to last
      else if (strcmp(chr2, "Y") == 0)
      {
        return true;
      }
      else if (strcmp(chr1, "Y") == 0)
      {
        return false;
      }
      // X just before Y
      else if (strcmp(chr2, "X") == 0)
      {
        return true;
      }
      else if (strcmp(chr1, "X") == 0)
      {
        return false;
      }

      // otherwise, in numerical order
      int numChr1 = toInt(chr1);
      int numChr2 = toInt(chr2);
      if (numChr1 < numChr2)
      {
        return true;
      }
      else if (numChr1 > numChr2)
      {
        return false;
      }
      else
      {
        return false;
      }
    }
    // chromosome names were same, so use position
    else if (pb1->chrBegin < pb2->chrBegin)
    {
      return true;
    }
    else
    {
      return false;
    }
  };
};

// This is for gene-oriented display
class GeneSummary
{

public:
  string name;
  vector<BlockSummary *> blocks; // vector of blocks overlapping this gene.
  vector<string> goTerms;
  bool isIgnored;
  GeneSummary(bool ignoreDefault) : isIgnored(ignoreDefault){};
};

// Genes comparator -- compares by best blocks in gene.
class GenesComparator
{
public:
  BlocksComparator bcomp;

  // constructor
  GenesComparator(bool isCat) : bcomp(BlocksComparator(isCat)){};

  bool operator()(const GeneSummary *pg1, const GeneSummary *pg2) const
  {
    BlockSummary *pb1 = pg1->blocks[0];
    BlockSummary *pb2 = pg2->blocks[0];
    return bcomp(pb1, pb2);
  };
};

// Globals
unordered_map<string, GeneSummary *> geneTable; // for gene-oriented interface

vector<BlockSummary *> blocks; // global vector of all blocks.

int traceFStat = false; // Debugging
// int traceFStat = true;

void filterGoTerms(char *fname, vector<string> terms)
{
  ColumnReader rdr(fname, (char *)"\t");
  int numtoks;
  vector<string>::iterator startT = terms.begin();
  vector<string>::iterator endT = terms.end();
  while ((numtoks = rdr.getLine()) >= 0)
  {
    string geneName = rdr.getToken(0);
    unordered_map<string, GeneSummary *>::iterator it = geneTable.find(geneName);
    if (it == geneTable.end())
    {
      continue;
    }
    for (int i = 1; i < numtoks; i++)
    {
      if (find(startT, endT, rdr.getToken(i)) != endT)
      {
        geneTable[geneName]->isIgnored = false;
      }
      geneTable[geneName]->goTerms.push_back(rdr.getToken(i));
    }
  }
}

// Some vector arithmetic.

// destroys first argument (like +=)
void addVectors(vector<float> &v1, vector<float> &v2)
{
  if (v1.size() != v2.size())
  {
    cout << "addVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << endl;
    exit(1);
  }
  vector<float>::iterator vend = v1.end();
  vector<float>::iterator vit2 = v2.begin();
  for (vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
  {
    *vit1 += *vit2;
    vit2++;
  }
}

void subtractVectors(vector<float> &v1, vector<float> &v2)
{
  if (v1.size() != v2.size())
  {
    cout << "subtractVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << endl;
    exit(1);
  }
  vector<float>::iterator vend = v1.end();
  vector<float>::iterator vit2 = v2.begin();
  for (vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
  {
    *vit1 -= *vit2;
    vit2++;
  }
}

//
float dotVectors(vector<float> &v1, vector<float> &v2)
{
  float result = 0.0;
  if (v1.size() != v2.size())
  {
    cout << "dotVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << endl;
    exit(1);
  }
  vector<float>::iterator vend = v1.end();
  vector<float>::iterator vit2 = v2.begin();
  for (vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
  {
    result += (*vit1) * (*vit2);
    vit2++;
  }
  return result;
}

// multiply by scalar.  Destroys first argument.
void scaleVector(vector<float> &v1, float c)
{
  vector<float>::iterator vend = v1.end();
  for (vector<float>::iterator vit = v1.begin(); vit < vend; vit++)
  {
    *vit *= c;
  }
}

// convert digits 0-9 in string to \000..\011 (but leave '?' printable).
void makeUnprintable(char *pattern)
{
  char *p = pattern;
  while (*p != 0)
  {
    if (*p != '?')
    {
      *p -= '0';
    }
    p++;
  }
}

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
        std::cout << "GCC: "<< __VERSION__ << std::endl;
        std::cout <<"ghmap version: "<<__HAPLOMAP_VER__<<std::endl;
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

// read in the block summary file.
// If geneName is non-null, only read the SNPs for that gene.
void readBlockSummary(char *fname, char *geneName, bool ignoreDefault)
{
  ColumnReader rdr(fname, (char *)"\t");

  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    // file has "Chromosome\tBlockNum\tBlockstartSNPnum\tSize\tChrbegin\tChrend\tPattern\tgene1\tcodon1\t...\n"

    // If geneName is defined, only process the lines containing that gene name in the gene name list,
    // unless geneName = "*", in which case it reads all blocks.

    if (!geneName || (find(rdr.begin() + 7, rdr.end(), geneName) != rdr.end()))
    {

      BlockSummary *pBlock = new BlockSummary(strdup(rdr.getToken(0).c_str()),
                                              toInt(rdr.getToken(1)),
                                              toInt(rdr.getToken(2)),
                                              toInt(rdr.getToken(3)),
                                              toInt(rdr.getToken(4)),
                                              toInt(rdr.getToken(5)),
                                              strdup(rdr.getToken(6).c_str()));
      // Add to blocks.
      blocks.push_back(pBlock);
      BlockSummary *pLastBlock = blocks.back();

      // For each gene, add gene summary to the gene table (if not there already),
      // and insert block in the block set.
      // This iterates over list of alternating gene name/codon tokens.
      for (vector<string>::iterator git = rdr.begin() + 7; git != rdr.end(); git += 2)
      {
        // *(git+1) is codon for gene name (*git).  Convert from string to bool using cond expr.
        pLastBlock->geneIsCodingMap[*git] = *(git + 1);
        unordered_map<string, GeneSummary *>::iterator entIt = geneTable.find(*git);

        if (entIt == geneTable.end())
        {
          // create new entry.
          // cout << "  Creating new gene summary for " << *git << endl;
          geneTable[*git] = new GeneSummary(ignoreDefault);
          geneTable[*git]->name = *git;
        }
        // Add block to GeneSummary
        //      cout << "  Adding block for " << *git << ": " << pLastBlock->blockIdx << endl;
        geneTable[*git]->blocks.push_back(pLastBlock);
      }

      // cout << "Gene names for block " << blocks.back().blockIdx << ": " << blocks.back().geneNames << endl;
    }
  }

  // Separate pass to fix up patterns.
  // patterns are all the same length.
  // FIXME: this causes a seg fault when there are no blocks (happens under funny conditions).
  numStrains = strlen(blocks[0]->pattern);

  for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
  {
    BlockSummary *block = blocks[blkIdx];
    makeUnprintable(block->pattern);
  }
}

// Order strain indices by decreasing phenotype value (or just group
// them, if categorical).
void sortStrainsByPheno(vector<vector<float>> &phenvec, vector<int> &strOrderVec)
{
  vector<int>::iterator stoEnd = strOrderVec.end();
  int i = 0;
  for (vector<int>::iterator stoIt = strOrderVec.begin(); stoIt != stoEnd; stoIt++)
  {
    *stoIt = i++;
  }
  // This will sort lexicographically, which is the right thing.
  // For categorical values, we just want equal values together.
  IndexComparator<vector<float>, less<vector<float>>> idxCompare(&phenvec);

  stable_sort(strOrderVec.begin(), strOrderVec.end(), idxCompare);
}

// Write summaries of the blocks.  The CGI script will read this and render it nicely in
// HTML.
void showGeneBlockSums(ostream &os, bool isCategorical, vector<BlockSummary *> &blocks, 
                       float cutoff, vector<int> &strOrderVec, vector<GeneSummary *> genesList)
{
  vector<string> genesOver; // all genes that have been already printed
  for (unsigned i = 0; i < blocks.size(); i++)
  {
    if (!isCutoff(isCategorical, cutoff, blocks[i]))
    {
      BlockSummary *pb = blocks[i];
      map<string, string> &geneIsCodingMap = pb->geneIsCodingMap;
      map<string, string>::iterator giit = geneIsCodingMap.begin();
      if (geneIsCodingMap.size() == 0)
      { // no genes in this block
        updateGeneIsInteresting(pb);
        os << "None\tNone\t";
        writeSortedPattern(os, pb->pattern, strOrderVec);
        os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) << "\t" << pb->effect << "\t";
        os << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t"; 
        os << pb->blockIdx << pb->blockStart << "\t" << pb->blockSize << endl;
      }
      else
      {
        for (; giit != geneIsCodingMap.end(); giit++)
        {
          string gname = giit->first;
          if (find(genesOver.begin(), genesOver.end(), gname) != genesOver.end())
          {
            continue;
          }
          genesOver.push_back(gname);
          for (vector<GeneSummary *>::iterator git = genesList.begin(); git != genesList.end(); git++)
          {
            if ((*git)->name == gname)
            {
              for (unsigned j = 0; j < ((*git)->blocks).size(); j++)
              {
                BlockSummary *pb_t = (*git)->blocks[j];
                updateGeneIsInteresting(pb_t);
                os << gname << "\t" << pb_t->geneIsInteresting[gname] << "\t";
                writeSortedPattern(os, pb_t->pattern, strOrderVec);
                os << "\t" << (isCategorical ? pb_t->FStat : pb_t->pvalue) << "\t" << pb_t->effect << "\t";
                os << pb_t->chrName << "\t" << pb_t->chrBegin << "\t" << pb_t->chrEnd << "\t";
                os << pb_t->blockIdx << pb_t->blockStart << "\t" << pb_t->blockSize;
                string upname = gname;
                upcase(upname);
                if (geneExprMap.find(upname) == geneExprMap.end())
                {
                  os << "\t-------------" << endl;
                }
                else
                {
                  os << "\t" << geneExprMap[upname] << endl;
                }
              }
              break;
            }
          }
        }
      }
    }
  }
}

// BY ADDITION
void writeGeneBlockSums(bool isCategorical, char *outputFileName, char *datasetName, vector<vector<float>> &phenvec, vector<BlockSummary *> &blocks, float pvalueCutoff)
{
  ofstream blockout(outputFileName);
  if (!blockout.is_open())
  {
    cout << "Open of file \"" << outputFileName << "\" failed: ";
    perror("");
    exit(1);
  }
  // Datasetname
  blockout << datasetName << endl;
  // sort strains by phenvec value
  vector<int> strOrderVec(numStrains); // will contain strain indices.
  sortStrainsByPheno(phenvec, strOrderVec);

  // output strain names.
  vector<int>::iterator stoEnd1 = strOrderVec.end();
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    blockout << strainAbbrevs.eltOf(str1);
    if (stoIt1 + 1 < stoEnd1)
    {
      blockout << "\t";
    }
  }
  blockout << endl;

  // output phenotype values.
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    blockout << phenvec[str1][0];
    if (stoIt1 + 1 < stoEnd1)
    {
      blockout << "\t";
    }
  }
  blockout << endl;

  GenesComparator gcomp(isCategorical);
  vector<GeneSummary *> genes;
  genes.reserve(geneTable.size());
  transform(geneTable.begin(), geneTable.end(), back_inserter(genes), std::bind(&std::unordered_map<string, GeneSummary *>::value_type::second, std::placeholders::_1 ));
            //select2nd<unordered_map<string, GeneSummary *>::value_type>());
           
  sort(genes.begin(), genes.end(), gcomp);
  showGeneBlockSums(blockout, isCategorical, blocks, pvalueCutoff, strOrderVec, genes);
}

void writeGeneBlockByBlocks(bool isCategorical, char *outputFileName, char *datasetName, vector<vector<float>> &phenvec, vector<BlockSummary *> &blocks, float pvalueCutoff)
{
  ofstream blockout(outputFileName);
  if (!blockout.is_open())
  {
    cout << "Open of file \"" << outputFileName << "\" failed: ";
    perror("");
    exit(1);
  }
  // Datasetname
  blockout << datasetName << endl;
  // sort strains by phenvec value
  vector<int> strOrderVec(numStrains); // will contain strain indices.
  sortStrainsByPheno(phenvec, strOrderVec);

  // output strain names.
  vector<int>::iterator stoEnd1 = strOrderVec.end();
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    blockout << strainAbbrevs.eltOf(str1);
    if (stoIt1 + 1 < stoEnd1)
    {
      blockout << "\t";
    }
  }
  blockout << endl;

  // output phenotype values.
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    blockout << phenvec[str1][0];
    if (stoIt1 + 1 < stoEnd1)
    {
      blockout << "\t";
    }
  }
  blockout << endl;

  showGeneBlockByBlocks(blockout, isCategorical, blocks, pvalueCutoff, strOrderVec);
}

void writeBlockSums(bool isCategorical, char *outputFileName,
                    char *datasetName, vector<vector<float>> &phenvec,
                    vector<BlockSummary *> &blocks, float pvalueCutoff)
{
  ofstream blockout(outputFileName);
  if (!blockout.is_open())
  {
    cout << "Open of file \"" << outputFileName << "\" failed: ";
    perror("");
    exit(1);
  }

  // Datasetname
  blockout << datasetName << endl;

  // sort strains by phenvec value

  vector<int> strOrderVec(numStrains); // will contain strain indices.
  sortStrainsByPheno(phenvec, strOrderVec);

  // output strain names.
  // FIXME:  Start using vecfuns written for Ravi microarray analysis
  vector<int>::iterator stoEnd1 = strOrderVec.end();
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    blockout << strainAbbrevs.eltOf(str1);
    if (stoIt1 + 1 < stoEnd1)
    {
      blockout << "\t";
    }
  }
  blockout << endl;

  // output phenotype values.
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    if (isCategorical)
    {
      blockout << catNames[str1];
    }
    else
    {
      blockout << phenvec[str1][0];
    }
    if (stoIt1 + 1 < stoEnd1)
    {
      blockout << "\t";
    }
  }
  blockout << endl;

  showBlockSums(blockout, isCategorical, blocks, pvalueCutoff, strOrderVec);
}

// Write gene-oriented summary.
void writeGeneSums(bool isCategorical, char *outputFileName,
                   char *datasetName, vector<vector<float>> &phenvec,
                   vector<BlockSummary *> &blocks, float cutoff, bool filterCoding)
{
  ofstream genesout(outputFileName);
  if (!genesout.is_open())
  {
    cout << "Open of file \"" << outputFileName << "\" failed: ";
    perror("");
    exit(1);
  }

  // Datasetname
  genesout << datasetName << endl;

  // sort strains by phenvec value

  vector<int> strOrderVec(numStrains); // will contain strain indices.
  sortStrainsByPheno(phenvec, strOrderVec);

  // output strain names.
  // FIXME:  Start using vecfuns written for Ravi microarray analysis
  //  ... or STL algorithms!
  vector<int>::iterator stoEnd1 = strOrderVec.end();
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    genesout << strainAbbrevs.eltOf(str1);
    if (stoIt1 + 1 < stoEnd1)
    {
      genesout << "\t";
    }
  }
  genesout << endl;

  // output phenotype values.
  for (vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
  {
    int str1 = *stoIt1;
    if (isCategorical)
    {
      genesout << catNames[str1];
    }
    else
    {
      genesout << phenvec[str1][0];
    }
    if (stoIt1 + 1 < stoEnd1)
    {
      genesout << "\t";
    }
  }
  genesout << endl;

  GenesComparator gcomp(isCategorical);

  // Copy genesTable values into a vector and sort using GenesComparator
  vector<GeneSummary *> genes;
  genes.reserve(geneTable.size());
  transform(geneTable.begin(), geneTable.end(), back_inserter(genes),
            std::bind(&std::unordered_map<string, GeneSummary *>::value_type::second, std::placeholders::_1 ));
            //select2nd<unordered_map<string, GeneSummary *>::value_type>());
  sort(genes.begin(), genes.end(), gcomp);

  // write them out.
  for (vector<GeneSummary *>::iterator git = genes.begin(); git != genes.end(); git++)
  {

    if ((*git)->isIgnored)
      continue; //simply don't write these genes
    BlockSummary *pBestBlock = (*git)->blocks[0];
    if (isCategorical)
    {
      if (pBestBlock->FStat < cutoff)
      {
        break;
      }
    }
    else
    {
      if (pBestBlock->pvalue > cutoff)
      {
        break;
      }
    }

    // write gene name and its coding bit.
    string &gname = (*git)->name;
    string ugname = gname; // gene names in expression data are upper case.
    upcase(ugname);

    // coding bit is subtle.  This iterates over all of the blocks that are < pvalue cutoff
    // (or > for FStat).  If any has a non-synonmous SNP, the bit is 1.  So, the coding
    // change may not be in the best block for the gene, and may change if the pvalue cutoff
    // is changed.
    //string isCodingStr = "";
    bool isCoding = false;
    bool hasInteresting = false;
    bool hasSpliceChange = false;
    //ofstream debug_log;
    //debug_log.open("debug.log",ios::app);
    for (vector<BlockSummary *>::iterator blit = (*git)->blocks.begin(); blit != (*git)->blocks.end(); blit++)
    {
      if (isCutoff(isCategorical, cutoff, *blit))
      {
        break;
      }
      else
      {
        //debug_log << gname << "\t" << (*blit)->geneIsCodingMap[gname] << endl;
        if ((*blit)->geneIsCodingMap[gname] != "0")
        { //if the gene had SNPs marked as NON_SYNONYMOUS_CODING, with <->, or as SPLICE_SITE, isCoding is true
          isCoding = true;
        }
        hasInteresting |= ((*blit)->numInteresting > 0);                                            //has a major amino acid change
        hasSpliceChange |= (((*blit)->geneIsCodingMap[gname]).find("SPLICE_SITE") != string::npos); //if "SPLICE_SITE" was in there that means that the gene had a splice change
      }
    }
    //New thing: -1 means not coding, 0 means coding but not important
    //Anything else is the number of important
    int codingCode = -1;
    if (isCoding)
      codingCode = 0;
    if (hasInteresting)
      codingCode = 1;
    if (hasSpliceChange)
      codingCode = 2; //if this is true, will overwrite codingCode=1
    updateGeneIsInteresting(pBestBlock);
    if ((isCoding || !filterCoding) || hasSpliceChange)
    {
      genesout << gname << "\t" << pBestBlock->geneIsInteresting[(*git)->name] << "\t";
      writeSortedPattern(genesout, pBestBlock->pattern, strOrderVec);
      genesout << "\t" << (isCategorical ? pBestBlock->FStat : pBestBlock->pvalue) << "\t" << pBestBlock->effect;
      genesout << "\t" << pBestBlock->chrName << "\t" << pBestBlock->chrBegin << "\t" << pBestBlock->chrEnd;

      // write gene expression values
      if (geneExprMap.find(ugname) == geneExprMap.end())
      {
        // no data for that gene name
        genesout << "\t-------------" << endl;
      }
      else
      {
        genesout << "\t" << geneExprMap[ugname] << endl;
      }
    }
  }
}

// Count the number of defined strains
int numDefinedStrains(char *pattern)
{
  int count = 0;
  for (int str1 = 0; str1 < numStrains; str1++)
  {
    if (pattern[str1] != '?')
    {
      count++;
    }
  }
  return count;
}

// returns maximum eq. class + 1.
int numHaplotypes(char *pattern)
{
  int numHap = -1;
  for (int str1 = 0; str1 < numStrains; str1++)
  {
    char hap = pattern[str1];
    if (pattern[str1] != '?' && numHap < hap)
    {
      numHap = hap;
    }
  }
  return numHap + 1;
}

/* Returns a score that represents the interestingness of codon changes*/
int scoreChanges(string str)
{
  int count = 0;
  size_t pos = 0;
  size_t endpos = 0;
  while (true)
  {
    pos = str.find("<", pos + 1);
    if (pos == string::npos)
      break;
    endpos = str.find(">", pos + 1);
    int aa1 = str[pos - 1] - 'A';
    int aa2 = str[endpos + 5] - 'A';
    if (AACLASSES[aa1] != AACLASSES[aa2])
    {
      count++;
      int X = 'X' - 'A';
      if (aa1 == X || aa2 == X)
        count += 2; // Higher weight if stop codon
    }
  }
  return count;
}

int interestingChanges(const map<string, string> &geneCodingMap)
{
  int changeCount = 0;
  for (map<string, string>::const_iterator git = geneCodingMap.begin();
       git != geneCodingMap.end(); git++)
  {
    if (git->second == "0" || git->second == "1")
      continue;
    changeCount += scoreChanges(git->second);
  }
  return changeCount;
}

// Mark each block as ignored unless it has a coding gene.
void filterCodingBlocks()
{
  for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
  {
    BlockSummary *pBlock = blocks[blkIdx];
    // look for a coding gene
    bool keep = false;
    map<string, string> &gicMap = pBlock->geneIsCodingMap;

    for (map<string, string>::iterator gicmit = gicMap.begin(); gicmit != gicMap.end(); gicmit++)
    {
      if ((*gicmit).second != "0")
      {
        keep = true;
        break;
      }
    }
    pBlock->isIgnored |= !keep; // |= as style to allow multiple filters.
  }
}

void filterEqualBlocks(vector<int> equalRegions)
{
  for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
  {
    BlockSummary *pBlock = blocks[blkIdx];
    // look for a coding gene
    bool keep = true;
    char allele = -1;
    for (unsigned int i = 0; i < equalRegions.size(); i++)
    {
      if (allele == -1)
      {
        allele = pBlock->pattern[equalRegions[i]];
        continue;
      }
      keep &= (allele == (pBlock->pattern)[equalRegions[i]]);
    }
    pBlock->isIgnored |= !keep; // |= as style to allow multiple filters.
  }
}

// Read a file of quantitative phenotypes.
void readQPhenotypes(char *fname, vector<vector<float>> &phenvec)
{
  ColumnReader rdr(fname, (char *)"\t");

  int numtoks;

  while ((numtoks = rdr.getLine()) >= 0)
  {
    // file has "Name\tAbbrev\n"
    if (numtoks != 2)
    {
      cout << "Warning: numtoks = " << numtoks << endl;
    }

    // FIXME: some unnecessary string copies
    string strain_abbrev = rdr.getToken(0);
    vector<float> qphen;
    qphen.push_back(toFloat(rdr.getToken(1)));
    //    int strIdx = strainAbbrevs.hasIndex(strain_abbrev);
    int strIdx = strainAbbrevs.addElementIfNew(strain_abbrev);
    if (strIdx < 0)
    {
      cout << "Undefined strain abbrev: " << strain_abbrev << endl;
    }
    phenvec[strIdx] = qphen;
  }
}

void setBlockStats()
{
  for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
  {
    BlockSummary *pBlock = blocks[blkIdx];
    //if(!pBlock->isIgnored) {
    pBlock->numHaplo = numHaplotypes(pBlock->pattern);
    pBlock->numInteresting = interestingChanges(pBlock->geneIsCodingMap);
  }
}
// Read a file of categorical phenotypes.
void readCPhenotypes(char *fname, vector<vector<float>> &phenvec)
{

  ColumnReader rdr(fname, (char *)"\t");
  Dynum<string> categories; // assigned the distinct categories consecutive indices, starting at 0.

  catNames.resize(numStrains);

  int numtoks;

  while ((numtoks = rdr.getLine()) >= 0)
  {
    // file has "Name\tAbbrev\n"
    if (numtoks != 2)
    {
      cout << "Warning: numtoks = " << numtoks << endl;
    }

    // FIXME: some unnecessary string copies
    string strain_abbrev = rdr.getToken(0);
    string catname = rdr.getToken(1);
    int strIdx = strainAbbrevs.addElementIfNew(strain_abbrev);
    catNames[strIdx] = catname;

    categories.addElementIfNew(catname);
  }

  // build category vectors.

  numCategories = categories.size();

  for (int strIdx = 0; strIdx < numStrains; strIdx++)
  {
    // build vector value for this category and store in phenvec.
    vector<float> cphen(categories.size(), 0.0F); // initialize to 0.
    cphen[categories.indexOf(catNames[strIdx])] = 1.0F;
    phenvec[strIdx] = cphen;
  }

  if (traceFStat)
  {
    cout << "Phenotype values" << endl;
    categories.dump();

    cout << "Phenotype vectors for strains" << endl;
    for (int str = 0; str < numStrains; str++)
    {
      cout << phenvec[str] << " ";
    }
    cout << endl;
  }
}

//read in equal class
void readEqualFile(char *fname, vector<int> &equalStrains)
{
  ColumnReader rdr(fname, (char *)"\t");
  //numStrains

  int numtoks;
  int i = 0;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    string eqclass = rdr.getToken(1);
    if (eqclass != "0")
      equalStrains.push_back(i);
    i++;
  }
}

//reads first column of tsv fname into vector vec
void readFileToVec(char *fname, vector<string> &vec)
{
  ColumnReader rdr(fname, (char *)"\t");
  //numStrains

  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    vec.push_back(rdr.getToken(0));
  }
}

// This version takes vector of vectors in phenotypes, so that it can handle normal
// and categorical data in the same code.
// Handling of categorical data is as described to me by Ming.
// Means are means of vectors.
// SSW is sum of squares of Euclidean distances of phenotype vectors to the phenotype mean.
// SSB is sum of squares of Euclidean distances of haplotype averages to global mean.
// This does not actually compute the p-value.  It stops with the F statistic.
void ANOVA(vector<vector<float>> &phenvec, char *pattern, float &FStat, float &pvalue, float &effect)
{
  int numHaplo = numHaplotypes(pattern);
  // array haplotype -> num strains in haplotype.
  vector<int> haploNum(numHaplo, 0);

  // array haplotype -> mean (vector<float>) for each haplotype
  vector<vector<float>> haploMean(numHaplo, vector<float>(numCategories, 0.0F));

  float numDefined = 0.0F;
  vector<float> sumDefined(numCategories, 0.0F);

  // Compute haplotype means
  for (int str1 = 0; str1 < numStrains; str1++)
  {
    int hap = pattern[str1];
    vector<float> &phen = phenvec[str1];
    if ('?' != hap)
    {
      numDefined++;
      addVectors(sumDefined, phen);
      haploNum[hap]++;

      addVectors(haploMean[hap], phen); // temporarily, the total, not mean.
    }
  }

  for (int hap = 0; hap < numHaplo; hap++)
  {
    scaleVector(haploMean[hap], 1.0 / haploNum[hap]); // get the mean
  }

  if (traceFStat)
  {
    cout << "numDefined = " << numDefined << ", sumDefined = " << sumDefined << endl;
    cout << "haploNum[] = [";
    for (int hap = 0; hap < numHaplo; hap++)
    {
      cout << haploNum[hap] << " ";
    }
    cout << "]" << endl;

    cout << "haploMean[] = [";
    for (int hap = 0; hap < numHaplo; hap++)
    {
      cout << haploMean[hap] << " ";
    }
    cout << "]" << endl;
  }

  float SSW = 0.0F;
  for (int str1 = 0; str1 < numStrains; str1++)
  {
    int hap = pattern[str1];
    if ('?' != hap)
    {
      vector<float> resid = phenvec[str1];
      subtractVectors(resid, haploMean[hap]);
      float tmpdot = dotVectors(resid, resid);

      SSW += tmpdot;
    }
  }

  // SSB -- between sum of squares (sum over haplotypes hapsize*(hapmean-mean)^2
  float SSB = 0.0;
  vector<float> &mean = sumDefined; // sumDefined will be the mean of all values.
  scaleVector(mean, 1.0 / numDefined);

  for (int hap = 0; hap < numHaplo; hap++)
  {
    vector<float> diff = haploMean[hap]; // copy so we don't destroy haploMeans
    subtractVectors(diff, mean);         // (haplotype mean) - mean
    float sq = haploNum[hap] * dotVectors(diff, diff);
    SSB += sq;
  }

  // my quick and dirty hack to penalize missing alleles.
  // simulates additional error for each missing value (but ignores
  // degrees of freedom).
  SSW += 1.1 * (numStrains - numDefined);

  // mean square within
  // df within is (numDefined-numHaplo)
  float dfW = numDefined - numHaplo; // degrees of freedom within
  float dfB = numHaplo - 1;          // degrees of freedom between.
  if (dfW == 0.0)
  {
    // This happens when numDefined = numHaplo, which occurs rarely when there are
    // lots of undefined strains and lots of haplotypes.
    // This causes errors in gsl, and I don't know what the right thing to do is,
    // so just punt.  We won't get a match for this.
    pvalue = 1.0;
    effect = 0.0;
    return;
  }
  float MSW = SSW / dfW;
  float MSB = SSB / dfB;

  // This formula is the same as Peltz.  So, I think I've seen two totally
  // different formulas for omega^2
  // WARNING: This divides by 0 if SSW is 0.
  // Which seems to work ok (F <- "inf").
  FStat = MSB / MSW;

  // out parameter for pvalue
  pvalue = (float)gsl_cdf_fdist_Q((double)FStat,
                                  (double)(numHaplo - 1),
                                  (double)(numDefined - numHaplo));

  // Genetic effect
  // Oh wow!  omega^2 is the "coefficient of determination"!
  // http://faculty.chass.ncsu.edu/garson/PA765/anova.htm#anova2
  effect = (float)((SSB - (numHaplo - 1) * MSW) / (SSW + SSB + MSW));

  if (traceFStat)
  {
    cout << "SSW = " << SSW
         << ", SSB = " << SSB
         << ", MSW = " << MSW
         << ", MSB = " << MSB
         << ", FStat = " << FStat
         << ", pval = " << pvalue
         << ", genetic effect = " << effect
         << endl;
  }

  // SST = SSW + SSB
  // (SSB/(k-1))/(SSW/(n-k)) -- k = numhaplotypes, n == numstrains (defined?)

  // Trying to correct between variations in notation/typos:  AP stats notes say:
  // F0 = MST(Between)/MSE(within)
  // MST(Between) = SST(between)/(df(Between)
  // MST(Within) = SST(Within)/(df(within))
  // Ok, so MSE(Between) = SSB/(k-1)
  // Ok, so MSE(Within) = SSW/(n-k)

  // "effect" "treatment" "between" seem to be roughly the same.
  // df(effect) is k-1
  // "error" "within" are equivalent.
  // df(error) is n-k

  // MSE - mean square error SSE/df(error)

  // omega^2 = (SSE k- df(effect)(MSE)) / (MSE + SST)
  // Peltz: omega^2 = (SSB - (k-1)*MSE)/(SST + MSE)
  //  INCONSISTENCY: I thought SSE was SSW
}

int main(int argc, char **argv)
{
  Options *opts = parseOptions(argc, argv);

  beginPhase("reading blocks summary file");
  readBlockSummary(opts->blocksFileName, opts->geneName, opts->goTermFile);
  endPhase();

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
  for (unordered_map<string, GeneSummary *>::iterator git = geneTable.begin(); git != geneTable.end(); git++)
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
