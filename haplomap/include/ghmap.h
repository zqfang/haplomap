#pragma once

#ifndef GHMAP_H
#define GHMAP_H

#include <cstring>
#include <string>
#include <unordered_map>
#include <functional>
#include "haploLib.h"

using namespace std;


int numHaplotypes(char *pattern); // forward decl.
int interestingChanges(const map<string, string> &geneCodingMap);


inline void upcase(string &str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

// Read the file of gene expression data
void readCompactGeneExpr(char *fname);

// renumber eqclasses in pattern so the increase from left to right
void writeSortedPattern(ostream &os, char *pattern, vector<int> &strOrderVec);


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
               int chrbeg, int chrend, char *pat);
  ~BlockSummary();

  string updateCodonScore(string str);
  friend void updateGeneIsInteresting(BlockSummary *pb);
  // print a line of the blocks file.
  // blockIdx	blockStart	blockSize	chromosome	begin	end	pattern	pval	effect	genename genehascoding ...
  friend void showBlockSum(ostream &os, bool isCategorical, BlockSummary *pb, vector<int> &strOrderVec);
  // BY addition, output as such:
  // gene	codon_flag	pattern	pval	effect	chromosome	begin	end	blockIdx	blockStart	blockSize	expression
  friend void showGeneBlockByBlock(ostream &os, bool isCategorical, BlockSummary *pb, vector<int> &strOrderVec);

  // Return true if pval is above cutoff or FStat is below it.
  friend bool isCutoff(bool isCategorical, float cutoff, BlockSummary *pBlock)
  {
    return (isCategorical) ? (pBlock->FStat < cutoff) : (pBlock->pvalue > cutoff);
  }
  friend void showBlockSums(ostream &os, bool isCategorical,
                            vector<BlockSummary *> &blocks, float cutoff, vector<int> &strOrderVec);
  friend void showGeneBlockByBlocks(ostream &os, bool isCategorical, vector<BlockSummary *> &blocks, float cutoff, vector<int> &strOrderVec);
};

inline void showIsCoding(map<string, string> geneIsCodingMap)
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
  bool operator()(const BlockSummary *pb1, const BlockSummary *pb2) const;
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


void filterGoTerms(char *fname, vector<string> terms);

// Some vector arithmetic.
// destroys first argument (like +=)
void addVectors(vector<float> &v1, vector<float> &v2);

void subtractVectors(vector<float> &v1, vector<float> &v2);

//
float dotVectors(vector<float> &v1, vector<float> &v2);

// multiply by scalar.  Destroys first argument.
void scaleVector(vector<float> &v1, float c);

// convert digits 0-9 in string to \000..\011 (but leave '?' printable).
void makeUnprintable(char *pattern);


// read in the block summary file.
// If geneName is non-null, only read the SNPs for that gene.
void readBlockSummary(char *fname, char *geneName, bool ignoreDefault);

// Order strain indices by decreasing phenotype value (or just group
// them, if categorical).
void sortStrainsByPheno(vector<vector<float>> &phenvec, vector<int> &strOrderVec);

// Write summaries of the blocks.  The CGI script will read this and render it nicely in
// HTML.
void showGeneBlockSums(ostream &os, bool isCategorical, vector<BlockSummary *> &blocks, 
                       float cutoff, vector<int> &strOrderVec, vector<GeneSummary *> genesList);

// BY ADDITION
void writeGeneBlockSums(bool isCategorical, char *outputFileName, char *datasetName,
        vector<vector<float>> &phenvec, vector<BlockSummary *> &blocks, float pvalueCutoff);

void writeGeneBlockByBlocks(bool isCategorical, char *outputFileName, char *datasetName,
        vector<vector<float>> &phenvec, vector<BlockSummary *> &blocks, float pvalueCutoff);

void writeBlockSums(bool isCategorical, char *outputFileName,
                    char *datasetName, vector<vector<float> > &phenvec,
                    vector<BlockSummary *> &blocks, float pvalueCutoff);

// Write gene-oriented summary.
void writeGeneSums(bool isCategorical, char *outputFileName,
                   char *datasetName, vector<vector<float> > &phenvec,
                   vector<BlockSummary *> &blocks, float cutoff, bool filterCoding);




/* Returns a score that represents the interestingness of codon changes*/
int scoreChanges(string str);

int interestingChanges(const map<string, string> &geneCodingMap);

// Mark each block as ignored unless it has a coding gene.
void filterCodingBlocks();

void filterEqualBlocks(vector<int> equalRegions);

// Read a file of quantitative phenotypes.
void readQPhenotypes(char *fname, vector<vector<float>> &phenvec);

void setBlockStats();
// Read a file of categorical phenotypes.
void readCPhenotypes(char *fname, vector<vector<float>> &phenvec);

//read in equal class
void readEqualFile(char *fname, vector<int> &equalStrains);

//reads first column of tsv fname into vector vec
void readFileToVec(char *fname, vector<string> &vec);

// This version takes vector of vectors in phenotypes, so that it can handle normal
// and categorical data in the same code.
// Handling of categorical data is as described to me by Ming.
// Means are means of vectors.
// SSW is sum of squares of Euclidean distances of phenotype vectors to the phenotype mean.
// SSB is sum of squares of Euclidean distances of haplotype averages to global mean.
// This does not actually compute the p-value.  It stops with the F statistic.
void ANOVA(vector<vector<float>> &phenvec, char *pattern, float &FStat, float &pvalue, float &effect);


// defined globals
/* use `extern` to declare global variable, to use it
 * You have to define these extern global variable in `main` file
 */
extern int numCategories; // default for non-categorical data.
extern const int AACLASSES[];
extern std::vector<string> catNames; // maps strIdx -> category name.
// Maps gene names to a string of A's, M's, and P's
extern std::unordered_map<string, string> geneExprMap;
// Globals
extern std::unordered_map<std::string, GeneSummary *> geneTable; // for gene-oriented interface
extern std::vector<BlockSummary *> blocks; // global vector of all blocks.
extern int traceFStat;

#endif