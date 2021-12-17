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
#include <unistd.h>
#include <cstring>
//for file lock to execute things in order
#include <fcntl.h>
#include <cerrno>
#include <cstdio>
#include <unordered_set>
#include <unordered_map>
#include "haplolib.h"

//using namespace std;



// Allocation of pattern strings.  all constructors should allocate.
// destructor should free never assign to the strings.  Do memcpy.
// Semantically, these should appear to be included in the structure,
// but we don't know the size to statically allocate it.

/*****************************************************************
 * Global variables and constants			  	 *
 *****************************************************************/
extern const int approxNumSNPs;
// If SNP is not defined for at least this many strains, ignore it.
extern int minDefined;
// FIXME: haploLimit should probably depend on number of strains
// tooLong should be configurable (options file?)
extern const int haploLimit;    // Maximum number of haplotypes to allow in a block.
extern const int tooLong; // A block cannot span more than 1 million base pairs

typedef std::vector<std::string> strvec; // vector of strings data type.

//// Map chromosome names to numerical indices.
//Dynum<string> chromosomes;
//// vector of all good SNPs in order of chromosome, position.
//vector<SNPInfo *> snpVec;

// Makes it easy to look up SNP info using SNP name.
//hash_map<string, SNPInfo *> snpMap(approxNumSNPs/2);
extern std::unordered_map<std::string, SNPInfo *> snpMap;
// Overlapping haplotype blocks, before we pick the best ones.
extern std::vector<HaploBlock> compoundHaploBlocks;
// set ordered by score.  Used as a priority queue.  Probably should
// use STL heap structure.
extern std::multiset<HaploBlock *, rCompareByScore> sortedCompoundHaploBlocks;
// vector of HaploBlocks chosen as best.
extern std::vector<HaploBlock *> chosenHaploBlocks;

/*****************************************************************
 * Functions and Methods				  	 *
 *****************************************************************/

// Read file of strains to use and index them.
int readDynumList(char *fname, Dynum<std::string> &dn);

// read file perlegen_b36_snp_vs_mmgene_091208.unl
// NEED THIS TO REGENERATE COMPACT FILE
// (Actually, not sure that the data is used.)
void readPerlegenSNPvsmgene(const char *fname);

// Read file perlegen_b36_strain_091208.unl
// Example: NES10630423_PWD/PhJ!PWD/PhJ!G!NES10630423!!
// Reading takes too long, so let's do this just once to generate a compact version.
// Maybe I should store all strains in this case.
void readPerlegenAlleleInfo(char *fname);

// Write in more compact format, basically an ascii dump of the snpMap.
// First line has a list of strains, in order
// Subsequent lines are
// SNPID\tchromosome\tposition\tallelestring\n
// allelestring is something like "AGGGA?AGAG", where alleles
// are in same order as strains on first line.
void writeAlleleInfoCompact(char *fname);

// read the compact format.
void readAlleleInfoCompact(char *fname);

// Read gene names and add to SNP info
void readSNPGeneNames(char *fname);

// A "good" SNP has no "D"'s, representing "gaps" (deletions?), has
// exactly two distinct alleles, and must have enough defined strains
// (minDefined elsewhere set to be >= 50% of alleles).
bool goodSNP(char *alleles, SNPInfo *SNPInfo);

// build the chromosome SNP lists.
void readChromosomeInfo(char *fname);

// read the Roche allele file and build associated data structures.
// Ignore irrelevant strains.
// Remove bad SNPs.
void readAlleleInfo(char *fname);

// set majorAllele to most frequent non-'?' char in alleles, minorAllele to second
// most frequent.  Code assumes there are only two alleles in string.  If a tie,
// pick the alphabetically first character as major allele.
void findMajorMinorAllele(char *alleles, char &majorAllele, char &minorAllele);

// Sort by chromosome first, then position
bool compareSNPs(SNPInfo *pSNP1, SNPInfo *pSNP2);

// Sort blocks by chromosome and then position
bool compareBlocksByPosition(HaploBlock *pHB1, HaploBlock *pHB2);
// Sort the SNP infos by chromosome, then by position.
void filterAndSortSNPs();

// Convert an allele string (e.g., "?ACCAT?") to a pattern (e.g., "?01102?").
// astr is allele string.  It is assumed not to have any 'D's at this point.
// Pattern is a pointer to a char array of length at least numStrains, which is
// written into.
// Returns true if pattern contains ?s
// FIXME: By convention, destination should be first.
// OBSOLETE: Returns number of haplotypes.
// FIXME: could uniquify pattern by hashing.  EASY!
int allelesToPattern(char *astr, char *pattern);

char minOfRowOrCol(char *table, int roworcol, int haploLimit, bool row);

int numHaplotypesFromAlleles(char *allelePattern);

inline double scoreBlock(int blocksize, int numHaplotypes)
{
  // (blocksize) / 2^(Number of haplotypes)
  return ((double)blocksize) / ((double)(1 << (numHaplotypes)));
}

int countHaplotypes(char *pattern);

// Given a pattern, renumber the colors so that they increase from left-to-right.
// This writes into it's argument.
void normalizePattern(char *pattern);

// Merge into "merge" array pattern eq classes of str1 in block of
// SNPs in snpVec starting at blockstart of size blocksize.
// Note: this uses same string encoding as pattern, but length is block size
// not numStrains (it's a column, not a row).  The purpose is to avoid
// combining columns that are individually compatible with the first
// column, but not all compatible with each other.
void mergePattern(char *merge, int blockstart, int blocksize, int str2);

// Given a block of SNPs within snpVec (blockstart is index of first SNP,
// blocksize is number of SNPs in block) and two strain indices, str1 and str2,
// check whether str1 and str2 have equal pattern chars, disregarding '?' entries.
// Merged is all combined columns in eq class so far.
bool strainsAreCompatible(char *merged, int blockstart, int blocksize, int str2);

bool qMarksInBlock(int blockstart, int blocksize);

int combinePatternNoQs(char *combined, int blockstart, int blocksize, int haploLimit);

// Greedy combine Patterns
int combinePatterns(char *combined, int blockstart, int blocksize, int haploLimit);

// Find highest-scoring block, starting at blockstart, only up to size maxSize
// and not exceeding haploLimit number of haplotypes.  Stores best block in by-ref first
// parameter.
void findBestBlock(HaploBlock &bestBlock, int blockstart, int haploLimit, int minSize, int maxSize);

// Search for haplotype increase after power-of-two block
// This is functionally equivalent to findBestBlock if applied after a max score power-of-two block.
// minSize is the power-of-two block size -- we search after this.
// haploLimit in this case is the number of haplotypes in the power-of-two block.
void ffbBinSearch(HaploBlock &bestBlock, int blockstart, int haploLimit, int minSize, int maxSize);

// Find highest-scoring block, starting at blockstart, only up to size maxSize
// and not exceeding haploLimit number of haplotypes.  Use size-doubling tricks
// to make it fast.
// WARNING: minSize is ignored.
void findBestBlockFast(HaploBlock &bestBlock, int blockstart, int haploLimit, int minSize, int maxSize );

// tests whether block is too big: too many base pairs, off end of
// chromosome.  Precondition: Doesn't go past end of snpVec.
bool blockIsTooBig(int blockstart, int blocksize, int chrIdx, int firstLoc);

// Find the longest block beginning at blockstart of 1<= minsize <= size <= blocksize that has
// haploSize number of haplotypes.  blocksize should be legitimate -- it should
// not cross too large a gap between SNPs or cross over a chromosome boundary.
// Allocates and returns the haplolock with the right pattern.
// If it returns true, pHaploBlock is properly filled in with pattern, size, etc.
// Otherwise, pHaploBlock is not written.
bool fmbBinSearch(HaploBlock *pHaploBlock, int blockstart, int minSize, int blocksize, int haploSize, int chrIdx, int firstLoc);

// From a given start SNP find the longest block for a specified
// number of haplotypes (haploSize), while crossing too many
// base pairs or chromosome boundary.  Return a pointer to the HaploBlock,
// or NULL if no such block exists (when there is a jump of more than 1 in
// the number of haplotypes).
// Method:  Double block sizes until haploLimit is exceeded.  Then do binary
// search for longest block with specified number of haplotypes.  There may
// be no such block (number of haplotypes may increase by more than one when a
// SNP is added).
HaploBlock *findMaximalBlock(int blockstart, int haploSize, int minSize);

// Find best block starting at each SNP.
void findCompoundBlocks(int haploLimit);

// For each haplotype size, find the maximal length blocks.
// *** Add trace functions.
void findAllMaximalBlocks(int haploLimit, int minSNPBlocks);

bool compareBySNP(HaploBlock *phb1, HaploBlock *phb2);

// Find best non-overlapping blocks.  Pick highest-score block, next highest-scoring
// that does not overlap, etc.
void chooseBlocks();

// print a block and the SNPs in it.
void showBlockSNPs(HaploBlock &hb);

void showBlocks(std::vector<HaploBlock> &hbv);

// Show a vector of pointers to blocks, instead of blocks. (Badly named.)
void showBlocksSNPs(std::vector<HaploBlock *> &hbs);

// strain display structure, used so we can put strains in a good
// order for HTML.
struct StrainDisp
{
  int strIdx;
  // strains in same haplotype.
  int numStrainsInSameHaplotype;
  // haplotype number (from pattern)
  int minStrainInSameHaplotype;
  // number of defined SNPs in block for this strain (non-'?' in column)
  int numSNPsDefined;
  // for debugging
  int haplotype;
  // constructor
  StrainDisp() : strIdx(0),
                 numStrainsInSameHaplotype(0),
                 minStrainInSameHaplotype(numStrains),
                 numSNPsDefined(0),
                 haplotype('?'){};

  friend std::ostream &operator<<(std::ostream &os, const StrainDisp &sd);
};

// comparison function to get in right order.
// haplotypes in decreasing order of number of strains.
// w/in haplotype, put most defined ones first?
//   haplotype -> straincount array.
//   strIdx (order in strain_index file) to break ties.
//  FIXME: Why can't I used references for sd1, sd2?
bool compareStrainDisps(StrainDisp sd1, StrainDisp sd2);

// Write out html for one block.
void writeBlock(std::ostream &hs, int blkIdx, int haploLimit);

// Writes HTML pages for a chromosome.  Since SNPs for multiple
// chromosomes may be in same vector and chosenBlocks may have the
// same property, this takes blkIdx, which is the index of the
// chosenBlock where the chromosome starts, and it returns the
// beginning of the next chromosome (or size+1 if at end). Returns
// first block of next chromosome (chosenHaploBlocks.size() when
// done).

int writeChromosome(const char *dirname, size_t blkIdx);

void writeHTML(char *dirName);

// write out list of gene names for a single block.
void writeBlockGeneNames(std::ofstream &os, HaploBlock *pHB);

// Write block summary for a single chromosome.
int writeChrBlockSummary(std::ofstream &os, size_t blkIdx, int minBlockSNPs);


void writeBlockSummary(char *fileName, int minBlockSNPs);

// Writes SNPs that were actually used in the block, in order of chromosome/position.
// This is specific to the set of strains and minBlockSNPs, because the good SNPs depend on the strains.
// This may include some SNPs in blocks that were discarded for being too small.
// FORMAT:  chrom	Pos	SNPID	Pattern		geneName1	codon1	geneName2 ...
void writeBlockSNPs(char *fname, int minBlockSNPs);

bool isCompatiblePattern(char *p1, char *p2);

// update stats after we have a chosenHaploBlocks vector (in case we do it after
// every chromosome).
void updateStats();

void reportStats();

