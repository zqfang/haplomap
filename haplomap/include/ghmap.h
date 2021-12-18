#pragma once

#ifndef GHMAP_H
#define GHMAP_H

#include <cstring>
#include <string>
#include <unordered_map>
#include <functional>
#include "haplolib.h"


// summary of a block, read for the file.
class BlockSummary
{
  public:
    char *chrName;
    int blockIdx;  // block id
    int blockStart; // snpvec index (0-based), output from eblocks -p
    int blockSize;  // snpvec size
    int chrBegin;
    int chrEnd;
    char *pattern;
    bool isIgnored; // for blocks that have been filtered out.
    float FStat;
    float pvalue;
    float FDR;
    float effect;
    float relFStat;
    float relPvalue;
    float relFDR;
    bool relReject;
    int numHaplo;
    int numInteresting;
    int numStrains;
    std::map<std::string, std::string> geneIsInteresting; // gene name -> codon change, by gene BY
    std::map<std::string, std::string> geneIsCodingMap;   // gene name -> coding strings, eg. TGA/X<->GGA/R
    // constructor
    BlockSummary(const char *chrnm, int num, int start, int size,
                int chrbeg, int chrend, const char *pat);
    ~BlockSummary();
    std::string updateCodonScore(std::string str);
    void updateGeneIsInteresting();
    void showIsCoding();
    // Return true if pval is above cutoff or FStat is below it.
    bool isCutoff(bool isCategorical, float cutoff);
    // print a line of the blocks file.
    // blockIdx	blockStart	blockSize	chromosome	begin	end	pattern	pval	effect	genename genehascoding ...
    friend void showBlockSum(std::ostream &os, bool isCategorical, BlockSummary *pb, std::vector<int> &strOrderVec);
    // BY addition, output as such:
    // gene	codon_flag	pattern	pval	effect	chromosome	begin	end	blockIdx	blockStart	blockSize	expression
    friend void showGeneBlockByBlock(std::ostream &os, bool isCategorical, BlockSummary *pb, std::vector<int> &strOrderVec);
    friend void showBlockSums(std::ostream &os, bool isCategorical,
                              std::vector<BlockSummary *> &blocks, float cutoff, std::vector<int> &strOrderVec);

private:
    void interestingChanges();
    int scoreChanges(std::string str);
    /// returns maximum eq. class + 1.
    int numHaplotypes();
    /// convert digits 0-9 in string to \000..\011 (but leave '?' printable).
    /// Mark: make ascii starts from 0, not '0'. Useful for downstream ANOVA analysis (a.k.a grouping number starts from 0)
    void makePatternUnprintable();
    /// return a printable Patten
    char * showPatternPrintable();
    /// whether a string is number
    bool isNumber(const std::string & str);
};


// This is to sort blocks in the best order.  First, sort by p-value
// or F statistic, depending on whether the phenotype is quantitative
// or categorical, then by chromosome, then by position.
class BlocksComparator
{
public:
  bool isCategorical;

  // constructor
  BlocksComparator(bool isCat);
  bool operator()(const BlockSummary *pb1, const BlockSummary *pb2) const;
};

// This is for gene-oriented display
class GeneSummary
{
public:
  std::string name;
  std::vector<BlockSummary *> blocks; // vector of blocks overlapping this gene.
  std::vector<std::string> goTerms;
  bool isIgnored;
  GeneSummary(bool ignoreDefault);
  ~GeneSummary();
};

// Genes comparator -- compares by best blocks in gene.
class GenesComparator
{
public:
  BlocksComparator bcomp;

  // constructor
  GenesComparator(bool isCat);
  bool operator()(const GeneSummary *pg1, const GeneSummary *pg2) const;
};



// upase strings
inline void upcase(std::string &str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

// Read the file of gene expression data
void readCompactGeneExpr(char *fname);

// 
void filterGoTerms(char *fname, std::vector<std::string> terms);

// read in the block summary file.
// If geneName is non-null, only read the SNPs for that gene.
void readBlockSummary(char *fname);

// Order strain indices by decreasing phenotype value (or just group
// them, if categorical).
// void sortStrainsByPheno(std::vector<std::vector<float>> &phenvec, std::vector<int> &strOrderVec);

// renumber eqclasses in pattern so the increase from left to right
void writeSortedPattern(std::ostream &os, char *pattern, std::vector<int> &strOrderVec);

// Write summaries of the blocks.  The CGI script will read this and render it nicely in HTML.
void showGeneBlockSums(std::ostream &os, bool isCategorical, std::vector<BlockSummary *> &blocks, 
                        float cutoff, std::vector<int> &strOrderVec, std::vector<GeneSummary *> genesList);
// BY ADDITION
void showGeneBlockByBlocks(std::ostream &os, bool isCategorical, std::vector<BlockSummary *> &blocks, 
                                    float cutoff, std::vector<int> &strOrderVec);

void writeGeneBlockSums(bool isCategorical, char *outputFileName, char *datasetName,
        std::vector<std::vector<float>> &phenvec, std::vector<BlockSummary *> &blocks, float pvalueCutoff);

void writeGeneBlockByBlocks(bool isCategorical, char *outputFileName, char *datasetName,
        std::vector<std::vector<float>> &phenvec, std::vector<BlockSummary *> &blocks, float pvalueCutoff);

void writeBlockSums(bool isCategorical, char *outputFileName, char *datasetName, 
                    std::vector<std::vector<float> > &phenvec,
                    std::vector<BlockSummary *> &blocks, float pvalueCutoff);

// Write gene-oriented summary.
void writeGeneSums(bool isCategorical, char *outputFileName,
                   char *datasetName, std::vector<std::vector<float> > &phenvec,
                   std::vector<BlockSummary *> &blocks, float cutoff, bool filterCoding);

// Write output headers
class GhmapWriter 
{
private:
    std::string _dataset_name;
    bool _isCategorical;

public:
    GhmapWriter(char *outputFileName, char *datasetName, bool isCategorical);
    ~GhmapWriter();
    void sortStrainsByPheno(std::vector<std::vector<float>> &phenvec, std::vector<int> &strOrderVec);
    void writeSortedPattern(char *pattern, std::vector<int> &strOrderVec); // NOTE: the input pattern has been makeUnprintable
    void writeStrainNameAndValue(std::vector<std::vector<float>> &phenvec, std::vector<int> &strOrderVec); // strainName and stainValues
    void writeExpressionNames(std::vector<std::string> &exprOrderVec);
    void writeHeaders(std::vector<std::string> &header);

    std::ofstream os;
};


// Mark each block as ignored unless it has a coding gene.
void filterCodingBlocks();

void filterEqualBlocks(std::vector<int> equalRegions);

// Read a file of quantitative phenotypes.
void readQPhenotypes(char *fname, std::vector<std::vector<float>> &phenvec, Dynum<std::string> & strainAbbrevs);

// Read a file of categorical phenotypes.
void readCPhenotypes(char *fname, std::vector<std::vector<float>> &phenvec, Dynum<std::string> & strainAbbrevs);

//read in equal class
void readEqualFile(char *fname, std::vector<int> &equalStrains);

//reads first column of tsv fname into vector vec
void readFileToVec(char *fname, std::vector<std::string> &vec);

// defined globals
/* use `extern` to declare global variable, to use it
 * You have to define these extern global variable in `main` file
 */
extern int numCategories; // default for non-categorical data.
extern const int AACLASSES[];
extern std::vector<std::string> catNames; // maps strIdx -> category name.
// Maps gene names to a string of A's, M's, and P's
extern std::unordered_map<std::string, std::string> geneExprMap;
// Gene expression header lines
extern std::vector<std::string> geneExprHeader;
// Globals
extern std::unordered_map<std::string, GeneSummary *> geneTable; // for gene-oriented interface
extern std::vector<BlockSummary *> blocks; // global vector of all blocks.
extern int traceFStat;

#endif