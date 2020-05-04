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
#include <unistd.h>
#include <cstring>
#include "haploLib.h"
//for file lock to execute things in order
#include <fcntl.h>
#include <cerrno>
#include <cstdio>
#include <unordered_set>
#include <unordered_map>

#ifndef __VERSION__
#def __VERSION__ 0.1
#endif

using namespace std;

// For setting initial datastructure sizes
const int approxNumSNPs = 8000000;

// debug/trace flags
//int  debugSNPMax = 1000;	// just do the first debugSNPMax SNPs. (-1 == all SNPs)
int debugSNPMax = -1;

bool traceFBB = false;
//bool traceFBB = true;
bool traceCombinePatterns = false;
//bool traceCombinePatterns = true;
bool traceChooseBlocks = false;
//bool traceChooseBlocks = true;
bool traceWriteHTML = false;

bool useFastFBB = true;
// bool useFastFBB = false;

// If SNP is not defined for at least this many strains, ignore it.
int minDefined;
int totalStrains = -1; // number of strains, not just relevant strains.

// Statistics counters
// 1) total number of blocks. 2) average number of SNPs per block, 3)
// average number of haplotypes per block, 4) number of small blocks (blocks
// that have only 1, 2 or 3 SNPs), and 5) number of SNPs in small blocks.

int totalBlocks = 0;
int totalSNPs = 0;
int totalHaploTypes = 0;
int numSmallBlocks = 0;
int totalSNPsInSmallBlocks = 0;

// Forward decls
int allelesToPattern(char *astr, char *pattern);

// Copy a pattern string.
char *dupPattern(char *pattern)
{
  char *newpat = (char *)malloc(numStrains);
  memcpy(newpat, pattern, numStrains);
  return newpat;
}

// Everything we need to know about a SNP.
class SNPInfo
{

  char *initAlleles()
  {
    if (numStrains < 0)
    {
      cout << "FATAL: Cannot allocate allele strings until strains_index.txt has been read." << endl;
      exit(1);
    }
    char *tmpAlleles = (char *)malloc(numStrains + 1);
    memset(tmpAlleles, '?', numStrains);
    tmpAlleles[numStrains] = '\000';
    return tmpAlleles;
  };

  // has to be called after initAlleles.
  char *initPattern()
  {
    char *tmpPattern = (char *)malloc(numStrains);
    memset(tmpPattern, '?', numStrains);
    return tmpPattern;
  };

public:
  friend ostream &operator<<(ostream &os, const SNPInfo &s);

  // default constructor
  SNPInfo() : alleles(initAlleles()), frozen(false), used(false), qMarks(true) { pattern = initPattern(); };

  // constructor for stuff in chr_info_perl file.
  //  SNPInfo(string n, int ci, int pos) : name(n), chrIdx(ci), position(pos) { alleles = initAlleles(); };
  SNPInfo(string n, int ci, int pos) : name(n), chrIdx(ci), position(pos),
                                       alleles(initAlleles()),
                                       frozen(false), used(false), qMarks(true)
  {
    pattern = initPattern();
  };

  // copy constructor
  SNPInfo(SNPInfo const &snpInfo) : name(snpInfo.name), chrIdx(snpInfo.chrIdx),
                                    position(snpInfo.position), alleles(strdup(snpInfo.alleles)),
                                    geneCodonMap(snpInfo.geneCodonMap),
                                    frozen(snpInfo.frozen), used(false), qMarks(true)
  {

    cerr << "Illegal call to SNPInfo copy constructor" << endl;
    exit(1);
    //     if (frozen) {
    //       cerr << "Copy constructor applied to frozen SNPInfo:" << this << endl;
    //     }
    //     cout << "Pattern copy: old " << hex << (size_t) snpInfo.pattern << ", new "
    // 	 << (size_t) pattern << endl;
    //     pattern = dupPattern(snpInfo.pattern);
  }

  // Overload assignment to memcpy, so we don't have to worry about string copying.
  SNPInfo &operator=(const SNPInfo &si)
  {
    cerr << "Illegal call to SNPInfo assignment operator." << endl;
    exit(1);
    if (this != &si)
    {
      if (frozen)
      {
        cerr << "Assignment of frozen SNPInfo:" << this << endl;
      }
      name = si.name;
      chrIdx = si.chrIdx;
      position = si.position;
      qMarks = si.qMarks;
      strcpy(alleles, si.alleles);
      // memcpy(pattern, si.pattern, numStrains);
    }
    return *this;
  };

  // destructor
  ~SNPInfo()
  {
    free(alleles);
    if (!frozen)
    {
      free(pattern);
    }
  }

  string name;
  int chrIdx;   // chromosome index
  int position; // location on chromosome
  // string of alleles "ACGTD?" "D" is deletion, "?" is unspecified.
  char *alleles;
  // string of haplotype numbers (not printable).
  char *pattern;
  map<string, string> geneCodonMap; // Map gene names to codons
  bool frozen;                      // when true, pattern is uniquified and should not be changed.
  bool used;                        // Is this SNP in a chosen block?
  bool qMarks;                      //Are there any unkowns in this SNP?

  char getAllele(int i) { return alleles[i]; }
  void setAllele(int i, char a) { alleles[i] = a; }
};

// Allocation of pattern strings.  all constructors should allocate.
// destructor should free never assign to the strings.  Do memcpy.
// Semantically, these should appear to be included in the structure,
// but we don't know the size to statically allocate it.

class HaploBlock
{
public:
  int start;   // Index in snpVec of first SNP in block.
  int size;    // number of good SNPs in block
  bool chosen; // Has it been chosen as a best compound block?

  // normalized allele pattern.  This uses '\0', '\1', etc.  and '?'
  // (obviously, assumes there are relatively few alleles).  NOTE:
  // These are not null-terminated strings, so strlen/strdup don't
  // work.
  char *pattern;
  double score; // block score

  void init()
  {
    pattern = (char *)malloc(numStrains);
    memset(pattern, '?', numStrains);
  };

  // default constructor
  HaploBlock() : chosen(false), score(0.0) { init(); };

  // copy constructor
  HaploBlock(HaploBlock const &hb) : start(hb.start), size(hb.size),
                                     chosen(hb.chosen), score(hb.score)
  {
    pattern = dupPattern(hb.pattern);
  };

  // Overload assignment to memcpy, so we don't have to worry about string copying.
  HaploBlock &operator=(const HaploBlock &hb)
  {
    if (this != &hb)
    {
      start = hb.start;
      size = hb.size;
      chosen = hb.chosen;
      score = hb.score;
      memcpy(pattern, hb.pattern, numStrains);
    }
    return *this;
  };

  // destructor
  ~HaploBlock() { free(pattern); };

  friend ostream &operator<<(ostream &os, const HaploBlock &s);
};

struct rCompareByScore
{
  bool operator()(const HaploBlock *phb1, const HaploBlock *phb2) const
  {
    return phb1->score > phb2->score;
  }
};

/*****************************************************************
 * Global variables and constants			  	 *
 *****************************************************************/

// FIXME: haploLimit should probably depend on number of strains
// tooLong should be configurable (options file?)

const int haploLimit = 5;    // Maximum number of haplotypes to allow in a block.
const int tooLong = 1000000; // A block cannot span more than 1 million base pairs

typedef vector<string> strvec; // vector of strings data type.

// Map chromosome names to numerical indices.
Dynum<string> chromosomes;

// Makes it easy to look up SNP info using SNP name.
//hash_map<string, SNPInfo *> snpMap(approxNumSNPs/2);
unordered_map<string, SNPInfo *> snpMap(approxNumSNPs / 2);
// vector of all good SNPs in order of chromosome, position.
vector<SNPInfo *> snpVec;

// Overlapping haplotype blocks, before we pick the best ones.
vector<HaploBlock> compoundHaploBlocks;

// set ordered by score.  Used as a priority queue.  Probably should
// use STL heap structure.
multiset<HaploBlock *, rCompareByScore> sortedCompoundHaploBlocks;

// vector of HaploBlocks chosen as best.
vector<HaploBlock *> chosenHaploBlocks;

// Test patterns for exact equality
inline bool eqPatterns(const char *pat1, const char *pat2)
{
  //  cout << "Pattern eq ";
  //  showPattern(pat1);
  //  cout << " = ";
  //  showPattern(pat2);
  bool result = memcmp(pat1, pat2, numStrains) == 0;
  //  cout << " result = " << result << endl;
  return result;
}

// hash function -- standard char* hash fails because of null termination
struct PatternHash
{
  size_t operator()(const char *s) const
  {
    //    cout << "Pattern hash: ";
    //    showPattern(s);
    // cout << " returns ";
    size_t h = 0;
    for (int i = 0; i < numStrains; i++)
    {
      h = 5 * h + s[i];
    }
    //    cout << h << endl;
    return h;
  }
};

// Hash table for pattern uniquification
struct PatternEq
{
  bool operator()(const char *pat1, const char *pat2) const
  {
    return eqPatterns(pat1, pat2);
  }
};

//typedef hash_set<char*, PatternHash, PatternEq> PatternSet;
typedef unordered_set<char *, PatternHash, PatternEq> PatternSet;
PatternSet patternUniqueTable;

/*****************************************************************
 * Functions and Methods				  	 *
 *****************************************************************/

// print a SNPInfo
ostream &operator<<(ostream &os, const SNPInfo &s)
{
  os << s.name << ", chromosome = " << chromosomes.eltOf(s.chrIdx)
     << ", location = " << s.position
     << ", alleles = " << s.alleles
     << ", pattern = ";
  showPattern(os, s.pattern);
  for (map<string, string>::const_iterator gmit = s.geneCodonMap.begin(); gmit != s.geneCodonMap.end(); gmit++)
  {
    os << ", " << (*gmit).first << " (" << (*gmit).second << ")";
  }
  os << endl;
  return os;
};

// print a HaploBlock
ostream &operator<<(ostream &os, const HaploBlock &hb)
{
  SNPInfo *pSNPInfo = snpVec[hb.start];
  os << "SNP " << hb.start << ": " << pSNPInfo->name
     << ", chromosome " << chromosomes.eltOf(pSNPInfo->chrIdx)
     << ", size " << hb.size
     << ", chosen  " << hb.chosen
     << ", pattern ";
  showPattern(hb.pattern);
  cout << ", score " << hb.score
       << endl;
  return os;
};

// for debugging
void showBlock(const HaploBlock &hb)
{
  cout << hb << endl;
}

void showBlocks(vector<HaploBlock *> &hbv)
{
  for (unsigned i = 0; i < hbv.size(); i++)
  {
    cout << *(hbv[i]) << endl;
  }
}

// Read file of strains to use and index them.
int readDynumList(char *fname, Dynum<string> &dn)
{
  ColumnReader rdr(fname, (char *)"\t");

  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    if (numtoks > 0)
    {
      dn.addElementIfNew(rdr.getToken(0));
    }
  }

  return dn.size();
}

// read file perlegen_b36_snp_vs_mmgene_091208.unl
// NEED THIS TO REGENERATE COMPACT FILE
// (Actually, not sure that the data is used.)
void readPerlegenSNPvsmgene(const char *fname)
{
  ColumnReader rdr(fname, (char *)"!");

  int numtoks;

  // example:	1!NES16340864!ENSMUSG00000073722!4931408C20Rik!!0!up  !96472!26728479!0!1!
  while ((numtoks = rdr.getLine()) >= 0)
  {

    string snpName = rdr.getToken(1);
    string geneName = rdr.getToken(3);
    int pos = atoi(rdr.getToken(8).c_str());
    string chr = rdr.getToken(10);

    int chrIdx = chromosomes.addElementIfNew(chr);

    SNPInfo *pSNPInfo = new SNPInfo(snpName, chrIdx, pos);
    if (!pSNPInfo)
    {
      cerr << "Fatal error: Could not allocation SNPInfo for " << snpName << "\n";
      abort();
    }

    // Insert a blank SNPEntry in the table with snpName
    // pair<hash_map<string, SNPInfo *>::iterator, bool> snpLookup =
    //   snpMap.insert(pair<string, SNPInfo *>(snpName, pSNPInfo));
    pair<unordered_map<string, SNPInfo *>::iterator, bool> snpLookup =
        snpMap.insert(pair<string, SNPInfo *>(snpName, pSNPInfo));
    // The file has multiple entries when there are multiple genes.  Just insert
    // the first, and don't check for duplicates.
    // if (!snpLookup.second) {
    //   cerr << "WARNING: Duplicate SNP names?" << localID << endl;
    // }
  }
}

// Read file perlegen_b36_strain_091208.unl
// Example: NES10630423_PWD/PhJ!PWD/PhJ!G!NES10630423!!
// Reading takes too long, so let's do this just once to generate a compact version.
// Maybe I should store all strains in this case.
void readPerlegenAlleleInfo(char *fname)
{
  ColumnReader rdr(fname, (char *)"!");

  int lineNum = 0;
  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {

    if (++lineNum % 10000 == 0)
    {
      cout << "." << flush;
    }

    string strain = rdr.getToken(1);
    string alleleStr = rdr.getToken(2);
    string SNPName = rdr.getToken(3);

    // hash_map<string, SNPInfo *>::iterator fit = snpMap.find(SNPName);
    unordered_map<string, SNPInfo *>::iterator fit = snpMap.find(SNPName);
    if (fit == snpMap.end())
    {
      // SNP name was not defined.  Warn and continue.
      //      cerr << "Warning: Undefined SNP name: " << SNPName << endl;
      continue;
    }

    SNPInfo *pSNPInfo = fit->second;

    int strIdx = strainAbbrevs.hasIndex(strain);

    // skip irrelevant strains.
    if (strIdx < 0)
    {
      continue;
    }

    char alleleChr = alleleStr[0];

    pSNPInfo->setAllele(strIdx, alleleChr);
  }
  cout << endl;
}

// Write in more compact format, basically an ascii dump of the snpMap.
// First line has a list of strains, in order
// Subsequent lines are
// SNPID\tchromosome\tposition\tallelestring\n
// allelestring is something like "AGGGA?AGAG", where alleles
// are in same order as strains on first line.
void writeAlleleInfoCompact(char *fname)
{
  ofstream cs(fname);
  if (!cs.is_open())
  {
    cerr << "Open of file \"" << fname << "\" failed: ";
    perror("");
    exit(1);
  }

  // write first line: strainAbbrevs (all strains, this time).
  for (int i = 0; i < numStrains; i++)
  {
    cs << strainAbbrevs.eltOf(i);
    if (i + 1 < numStrains)
    {
      cs << "\t";
    }
  }
  cs << endl;
  // iterate over hash table writing the rest of the strains.
  // not "pattern" and "used", though.
  // hash_map<string, SNPInfo *>::iterator send = snpMap.end();
  unordered_map<string, SNPInfo *>::iterator send = snpMap.end();
  // for (hash_map<string, SNPInfo *>::iterator sit = snpMap.begin(); sit != send; sit++) {
  //   SNPInfo * pSNPInfo = sit->second;
  //   cs << pSNPInfo->name << "\t" << chromosomes.eltOf(pSNPInfo->chrIdx) << "\t"
  //      << pSNPInfo->position << "\t" << pSNPInfo->alleles << endl;
  // }
  for (unordered_map<string, SNPInfo *>::iterator sit = snpMap.begin(); sit != send; sit++)
  {
    SNPInfo *pSNPInfo = sit->second;
    cs << pSNPInfo->name << "\t" << chromosomes.eltOf(pSNPInfo->chrIdx) << "\t"
       << pSNPInfo->position << "\t" << pSNPInfo->alleles << endl;
  }
}

// read the compact format.
void readAlleleInfoCompact(char *fname)
{
  ColumnReader rdr(fname, (char *)"\t");

  Dynum<string> allStrains;

  int numtoks;
  // get line with strains, and build list of all strains (relevant
  // strains depends on phenotypes, while these are all the good sequenced strains).
  numtoks = rdr.getLine();
  for (int aStrIdx = 0; aStrIdx < numtoks; aStrIdx++)
  {
    string strain = rdr.getToken(aStrIdx);
    allStrains.addElementIfNew(strain);
  }

  // read the lines with the SNP info
  while ((numtoks = rdr.getLine()) >= 0)
  {

    string snpName = rdr.getToken(0);
    // HACK:  if snp name ends with "_incon", remove that suffix before
    // looking it up in the gene name mapping.  This was Ming's way of annotating
    // SNPs that were inconsistent between the Sanger and Perlegen SNPs.
    // DILL 5/22/10 disabled this, because the full name is used to track gene names.
    if (false && snpName.compare(snpName.size() - 6, 6, "_incon") == 0)
    {
      //      cout << "removing _incon from: " << snpName << " -> ";
      snpName = snpName.substr(0, snpName.size() - 6);
      //      cout << snpName << endl;
    }

    string chr = rdr.getToken(1);
    int chrIdx = chromosomes.addElementIfNew(chr);
    int pos = atoi(rdr.getToken(2).c_str());
    string alleleStr = rdr.getToken(3);

    if (totalStrains < 0)
    {
      totalStrains = alleleStr.size();
    }

    // Only process SNP if it has at least ceil(totalStrains/2) defined.
    // We do this here, not in filterAndSortSNPs, because we have the total number
    // of strains and number defined among those here.
    int numDefined = totalStrains;
    string::iterator end = alleleStr.end();
    for (string::iterator it = alleleStr.begin(); it < end; it++)
    {
      if ('?' == *it)
      {
        numDefined--; //remove strains if SNP is '?'
      }
    }

    if (numDefined >= (totalStrains + 1) / 2)
    {
      SNPInfo *pSNPInfo = new SNPInfo(snpName, chrIdx, pos);
      if (!pSNPInfo)
      {
        cerr << "Fatal error: Could not allocation SNPInfo for " << snpName << "\n";
        abort();
      }

      // Insert a blank SNPEntry in the table with snpName
      //     pair<hash_map<string, SNPInfo *>::iterator, bool> snpLookup =
      // snpMap.insert(pair<string, SNPInfo *>(snpName, pSNPInfo));
      pair<unordered_map<string, SNPInfo *>::iterator, bool> snpLookup =
          snpMap.insert(pair<string, SNPInfo *>(snpName, pSNPInfo));

      if (!snpLookup.second)
      {
        cout << "WARNING: SNP " << snpName << " is duplicated in " << fname << endl;
        continue; // just skip duplicate one, although this should never happen.
      }

      // build up the allele string
      for (int strIdx = 0; strIdx < numStrains; strIdx++)
      {
        //string strain = relevantStrains.eltOf(strIdx);
        string strain = strainAbbrevs.eltOf(strIdx);
        // will report error if relevant strain was not in data.
        int aStrIdx = allStrains.indexOf(strain);
        char alleleChr = alleleStr[aStrIdx];
        pSNPInfo->setAllele(strIdx, alleleChr);
      }
      //cout << pSNPInfo->alleles << "\t" << snpName << endl;
    }
  }
}

// Read gene names and add to SNP info
void readSNPGeneNames(char *fname)
{
  ColumnReader rdr(fname, (char *)"\t");
  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    string snpName = rdr.getToken(0);
    string geneName = rdr.getToken(1);

    // cout << " geneName = " << geneName << endl;

    // hash_map<string, SNPInfo *>::iterator fit = snpMap.find(snpName);
    unordered_map<string, SNPInfo *>::iterator fit = snpMap.find(snpName);
    if (fit == snpMap.end())
    {
      // SNP name was not defined.  Warn and continue.
      //      cerr << "Warning: Undefined SNP name: " << snpName << endl;
      continue;
    }

    SNPInfo *pSNPInfo = fit->second;

    // rest of line is alternating gene names/codons

    for (vector<string>::iterator rit = rdr.begin() + 1; rit != rdr.end(); rit += 2)
    {
      // Store "<" in the map even if there is already something there.
      if (pSNPInfo->geneCodonMap.find(*rit) == pSNPInfo->geneCodonMap.end() || (*(rit + 1)).find("<") != string::npos)
      {
        pSNPInfo->geneCodonMap[*rit] = *(rit + 1);
      }
    }
  }
}

// A "good" SNP has no "D"'s, representing "gaps" (deletions?), has
// exactly two distinct alleles, and must have enough defined strains
// (minDefined elsewhere set to be >= 50% of alleles).
bool goodSNP(char *alleles, SNPInfo *SNPInfo)
{
  char firstAllele = '?';
  char secondAllele = '?';

  bool polymorphic = false;
  int numDefined = 0;

  for (int i = 0; i < numStrains; i++)
  {
    char allele = alleles[i];
    if (allele == 'D')
    {
      // cout << "  Eliminating SNP because of gap:\t" << snpInfo << endl;
      return false; // has a gap, exclude it.
    }
    if ('?' != allele)
    {
      // allele is defined.
      numDefined++;

      if ('?' == firstAllele)
      {
        // allele is first defined one.
        firstAllele = allele;
      }
      else if (allele != firstAllele)
      {
        // This defined allele is different from first one.
        polymorphic = true;
        if ('?' == secondAllele)
        {
          // just two alleles.
          secondAllele = allele;
        }
        else if (allele != secondAllele)
        {
          // cout << "Eliminating SNP because more than two alleles:\t" << snpInfo << endl;
          return false;
        }
        // else allele defined, firstAllele & secondAllele defined, allele == secondAllele
      }
      // else allele defined, firstAllele defined, allele == firstAllele
      // do nothing
    }
  }

  if (!polymorphic)
  {
    // only saw one allele
    // cout << "  Eliminating SNP  because not polymorphic:\t"  << snpInfo << endl;
    return false;
  }
  else if (numDefined >= minDefined)
  {
    return true;
  }
  else
  {
    // cout << "  Eliminating SNP  because of too few defined strains:\t" << snpInfo << endl;
    return false;
  }
}

// build the chromosome SNP lists.
void readChromosomeInfo(char *fname)
{
  ColumnReader rdr(fname, (char *)"\t");

  // parse chr_info_perl.txt
  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    // file has "SNPname\t\tchr\t\tpos\n"

    // FIXME: some unnecessary string copies
    string snpName = rdr.getToken(0);
    string chr = rdr.getToken(2);
    int chrIdx = chromosomes.indexOf(chr);
    int pos = toInt(rdr.getToken(4));

    // Insert a blank SNPEntry in the table with snpName
    SNPInfo *pSNPInfo = new SNPInfo(snpName, chrIdx, pos);
    if (!pSNPInfo)
    {
      cerr << "Fatal error: Could not allocation SNPInfo for " << snpName << "\n";
      abort();
    }

    // pair<hash_map<string, SNPInfo *>::iterator, bool> snpLookup =
    //   snpMap.insert(pair<string, SNPInfo*>(snpName, pSNPInfo));
    pair<unordered_map<string, SNPInfo *>::iterator, bool> snpLookup =
        snpMap.insert(pair<string, SNPInfo *>(snpName, pSNPInfo));
    if (!snpLookup.second)
    {
      cerr << "WARNING: Duplicate SNP names?" << snpName << endl;
    }
  }
}

// read the Roche allele file and build associated data structures.
// Ignore irrelevant strains.
// Remove bad SNPs.
void readAlleleInfo(char *fname)
{
  ColumnReader rdr(fname, (char *)"\t");
  // this has the weird double-tabs, too.

  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    if (numtoks != 5)
    {
      cerr << "Warning: numtoks = " << numtoks << endl;
    }

    string snpName = rdr.getToken(0);

    SNPInfo *pSNPInfo = snpMap[snpName];

    string strain = rdr.getToken(2);

    int strIdx = strainAbbrevs.hasIndex(strain);

    // skip irrelevant strains.
    if (strIdx < 0)
    {
      continue;
    }

    string alleleStr = rdr.getToken(4);
    char alleleChr = (alleleStr == "GAP") ? 'D' : alleleStr[0];

    pSNPInfo->setAllele(strIdx, alleleChr);
  }
}

// set majorAllele to most frequent non-'?' char in alleles, minorAllele to second
// most frequent.  Code assumes there are only two alleles in string.  If a tie,
// pick the alphabetically first character as major allele.
void findMajorMinorAllele(char *alleles, char &majorAllele, char &minorAllele)
{

  majorAllele = minorAllele = 0;

  int ACount = 0, CCount = 0, GCount = 0, TCount = 0;
  for (int strIdx = 0; strIdx < numStrains; strIdx++)
  {
    switch (alleles[strIdx])
    {
    case 'A':
      ACount++;
      break;
    case 'C':
      CCount++;
      break;
    case 'G':
      GCount++;
      break;
    case 'T':
      TCount++;
      break;
    }
  }

  int tot = ACount + CCount + GCount + TCount;
  if (ACount >= tot - ACount)
  {
    majorAllele = 'A';
  }
  else if (ACount > 0)
  {
    minorAllele = 'A';
  }

  if (CCount >= tot - CCount && majorAllele == 0)
  {
    majorAllele = 'C';
  }
  else if (CCount > 0)
  {
    minorAllele = 'C';
  }

  if (GCount >= tot - GCount && majorAllele == 0)
  {
    majorAllele = 'G';
  }
  else if (GCount > 0)
  {
    minorAllele = 'G';
  }

  if (TCount > tot - TCount && majorAllele == 0)
  {
    majorAllele = 'T';
  }
  else if (TCount > 0)
  {
    minorAllele = 'T';
  }
}

// Sort by chromosome first, then position
bool compareSNPs(SNPInfo *pSNP1, SNPInfo *pSNP2)
{
  if (pSNP1->chrIdx == pSNP2->chrIdx)
  {
    return (pSNP1->position < pSNP2->position);
  }
  else
  {
    return pSNP1->chrIdx < pSNP2->chrIdx;
  }
}

// Sort blocks by chromosome and then position
bool compareBlocksByPosition(HaploBlock *pHB1, HaploBlock *pHB2)
{
  SNPInfo *pSNP1 = snpVec[pHB1->start];
  SNPInfo *pSNP2 = snpVec[pHB2->start];
  return compareSNPs(pSNP1, pSNP2);
}

// Sort the SNP infos by chromosome, then by position.
void filterAndSortSNPs()
{
  // FIXME:  Be careful about goodSNP.  I don't think we want to just delete SNPs without
  // minStrains strains, because a SNP with 2 strains might break up a haplotype block --
  // and we might not want to miss that (unless few strains means the SNP is not reliable).
  // FIXME:  can we use a reference or pointer to reduce copies?
  int snpCount = 0;
  int goodSNPCount = 0;
  int badSNPCount = 0;

  ++snpCount;
  // FIXME: delete snpInfo's as we go through this?
  // for (hash_map<string, SNPInfo *>::iterator it = snpMap.begin(); it != snpMap.end(); it++)
  for (unordered_map<string, SNPInfo *>::iterator it = snpMap.begin(); it != snpMap.end(); it++)
  {
    SNPInfo *pSNPInfo = (*it).second;
    if (goodSNP(pSNPInfo->alleles, pSNPInfo))
    {
      // Initialize pattern (had to wait until all alleles were read)
      bool qMarks = allelesToPattern(pSNPInfo->alleles, pSNPInfo->pattern);
      pSNPInfo->qMarks = qMarks;
      pair<PatternSet::iterator, bool> p = patternUniqueTable.insert(pSNPInfo->pattern);
      if (!p.second)
      {
        // it was already there, so we didn't add it.
        free(pSNPInfo->pattern);
        pSNPInfo->pattern = *(p.first);
      }
      pSNPInfo->frozen = true;
      //if (!qMarks) // debug here
      snpVec.push_back(pSNPInfo);
      goodSNPCount++;
    }
    else
    {
      badSNPCount++;
    }
  }

  //  cout << "SNP filtering: good = " << goodSNPCount
  //       << ", bad = " << badSNPCount
  //       << ", total = " << snpCount << endl;

  // sort them with funky comparison function.
  sort(snpVec.begin(), snpVec.end(), compareSNPs);
  // just take the first debugSNPMax SNPs for debugging
  if (debugSNPMax > 0 && debugSNPMax < (int)snpVec.size())
  {
    snpVec.resize(debugSNPMax);
  }
}

// Convert an allele string (e.g., "?ACCAT?") to a pattern (e.g., "?01102?").
// astr is allele string.  It is assumed not to have any 'D's at this point.
// Pattern is a pointer to a char array of length at least numStrains, which is
// written into.
// Returns true if pattern contains ?s
// FIXME: By convention, destination should be first.
// OBSOLETE: Returns number of haplotypes.
// FIXME: could uniquify pattern by hashing.  EASY!
int allelesToPattern(char *astr, char *pattern)
{
  bool qMarks = false;
  int eqclass = 0;

  char map[128];
  map[(int)'A'] = '?';
  map[(int)'C'] = '?';
  map[(int)'G'] = '?';
  map[(int)'T'] = '?';
  map[(int)'?'] = '?';

  for (int strIdx = 0; strIdx < numStrains; strIdx++)
  {
    int a = astr[strIdx];
    if (a == '?')
    {
      qMarks = true;
      // no change to map.  Pattern[strIdx] stays '?'
    }
    else if (map[a] == '?')
    {
      // a not yet mapped.  Bind it.
      map[a] = eqclass;
      pattern[strIdx] = eqclass++;
    }
    else
    {
      pattern[strIdx] = map[a];
    }
  }
  return qMarks;
}

char minOfRowOrCol(char *table, int roworcol, int haploLimit, bool row)
{
  char result = 127;
  char entry;
  for (int i = 0; i < haploLimit; i++)
  {
    if (row)
    {
      entry = table[roworcol * haploLimit + i];
    }
    else
    {
      entry = table[i * haploLimit + roworcol];
    }
    if (entry < result)
    {
      result = entry;
    }
  }
  if (result < 0 || (result != '?' && result > haploLimit))
  {
    cerr << "BUG!!!" << endl;
  }
  return result;
}

int numHaplotypesFromAlleles(char *allelePattern)
{
  int Acount = 0;
  int Ccount = 0;
  int Gcount = 0;
  int Tcount = 0;
  for (int i = 0; i < numStrains; i++)
  {
    switch (allelePattern[i])
    {
    case 'A':
    {
      Acount++;
      break;
    }
    case 'C':
    {
      Ccount++;
      break;
    }
    case 'G':
    {
      Gcount++;
      break;
    }
    case 'T':
    {
      Tcount++;
      break;
    }
    }
  }
  return (Acount > 0) + (Ccount > 0) + (Gcount > 0) + (Tcount > 0);
}

double scoreBlock(int blocksize, int numHaplotypes)
{
  // (blocksize) / 2^(Number of haplotypes)

  return ((double)blocksize) / ((double)(1 << (numHaplotypes)));
}

int countHaplotypes(char *pattern)
{
  int numHaplo = 0;
  for (int i = 0; i < numStrains; i++)
  {
    char eqclass = pattern[i];
    if ('?' != eqclass)
    {
      numHaplo = (eqclass > numHaplo) ? eqclass : numHaplo;
    }
  }
  return numHaplo + 1;
}

// Given a pattern, renumber the colors so that they increase from left-to-right.
// This writes into it's argument.
void normalizePattern(char *pattern)
{
  // map from old to new.
  char cmap[haploLimit + 1];
  memset(&cmap, '?', haploLimit + 1);

  char eqclass = 0;

  for (int str1 = 0; str1 < numStrains; str1++)
  {
    // since this overwrites argument, be not to read the over-written value!
    int oldeq = pattern[str1];
    if (oldeq != '?')
    {
      if (cmap[oldeq] == '?')
      {
        // normalized equivalence class has not been assigned.
        cmap[oldeq] = eqclass++;
      }
      pattern[str1] = cmap[oldeq];
    }
  }
}

// Merge into "merge" array pattern eq classes of str1 in block of
// SNPs in snpVec starting at blockstart of size blocksize.
// Note: this uses same string encoding as pattern, but length is block size
// not numStrains (it's a column, not a row).  The purpose is to avoid
// combining columns that are individually compatible with the first
// column, but not all compatible with each other.
void mergePattern(char *merge, int blockstart, int blocksize, int str2)
{
  vector<SNPInfo *>::iterator blockBeginIt = snpVec.begin() + blockstart;
  vector<SNPInfo *>::iterator blockEndIt = blockBeginIt + blocksize;

  for (vector<SNPInfo *>::iterator snpIt = blockBeginIt; snpIt != blockEndIt; snpIt++)
  {
    int snpOffset = (snpIt - blockBeginIt); // position of current SNP w/in block.
    char *pat = (*snpIt)->pattern;
    char &chr1 = merge[snpOffset];
    char chr2 = pat[str2];
    if ('?' == chr1)
    {
      chr1 = chr2; // updates merge (chr1 is ref), if merge == ?, set to str2 pattern
    }
  }
}

// Given a block of SNPs within snpVec (blockstart is index of first SNP,
// blocksize is number of SNPs in block) and two strain indices, str1 and str2,
// check whether str1 and str2 have equal pattern chars, disregarding '?' entries.
// Merged is all combined columns in eq class so far.
bool strainsAreCompatible(char *merged, int blockstart, int blocksize, int str2)
{
  vector<SNPInfo *>::iterator blockBeginIt = snpVec.begin() + blockstart;
  vector<SNPInfo *>::iterator blockEndIt = blockBeginIt + blocksize;
  char *prevpat = NULL;
  for (vector<SNPInfo *>::iterator snpIt = blockBeginIt; snpIt != blockEndIt; snpIt++)
  {
    // FIXME: try fast exit when SNPs are identical
    char *pat = (*snpIt)->pattern;
    if (pat == prevpat)
    {
      continue; // pattern won't change anything.
    }
    int snpOffset = (snpIt - blockBeginIt); // position of current SNP w/in block.
    char chr1 = merged[snpOffset];
    char chr2 = pat[str2];
    if (chr1 != chr2 && '?' != chr1 && '?' != chr2)
    {
      return false;
    }
    prevpat = pat;
  }
  return true;
}

bool qMarksInBlock(int blockstart, int blocksize)
{
  vector<SNPInfo *>::iterator blockBeginIt = snpVec.begin() + blockstart;
  vector<SNPInfo *>::iterator blockEndIt = blockBeginIt + blocksize;
  for (vector<SNPInfo *>::iterator snpIt = blockBeginIt; snpIt != blockEndIt; snpIt++)
  {
    if ((*snpIt)->qMarks)
      return true;
  }
  return false;
}

int combinePatternNoQs(char *combined, int blockstart, int blocksize, int haploLimit)
{
  vector<SNPInfo *>::iterator blockBeginIt = snpVec.begin() + blockstart;
  vector<SNPInfo *>::iterator blockEndIt = blockBeginIt + blocksize;

  int firstInClassIndex = 0;
  int firstInClass = -1;
  int numAssigned = 0;
  for (int curClass = 0; curClass < haploLimit; curClass++)
  {
    for (int i = 0; i < numStrains; i++)
    {
      if (combined[i] != '?')
        continue; //because this has been assigned already;
      if (firstInClass == -1)
      {
        firstInClass = i;
        combined[i] = curClass;
        numAssigned++;
        continue;
      }
      combined[i] = curClass;
      numAssigned++;
      for (vector<SNPInfo *>::iterator snpIt = blockBeginIt; snpIt != blockEndIt; snpIt++)
      {
        int index = snpIt - blockBeginIt;
        char *pat = (*snpIt)->pattern;
        if (pat[i] != pat[firstInClass]) // each one compared to the first pattern in blocks
        {
          combined[i] = '?';
          numAssigned--;
          break;
        }
      }
    }
    if (firstInClass == -1)
    {
      //showPattern(combined);
      //cout << "from" << endl;
      //cout << endl;
      return curClass;
    }
    firstInClass = -1;
  }
  //cout << endl;
  if (numAssigned == numStrains)
    return haploLimit;
  return 0;
}

// Greedy combine Patterns
int combinePatterns(char *combined, int blockstart, int blocksize, int haploLimit)
{
  if (blocksize == 1)
  {
    memcpy(combined, snpVec[blockstart]->pattern, numStrains);
    return 2; //This is hacky: If there's only 1 or >2, this function will never be called
              /*int max = 0;
    for( i = 0; i < numStrains; i++) {
      char eqclass = combined[i];
      if(eqclass > max && eqclass!='?')
	max = eqclass;
    }
    return max + 1; */
  }
  if (traceCombinePatterns)
  {
    cout << "Combining patterns for " << blocksize << " SNPs starting at " << blockstart << endl;
  }
  if (traceCombinePatterns)
  {
    for (int snpIdx = blockstart; snpIdx < blockstart + blocksize; snpIdx++)
    {
      cout << "    ";
      showPattern(snpVec[snpIdx]->pattern);
      cout << endl;
    }
  }

  // Build vector of counts of defined entries in each column.
  // Flag columns consisting of all '?' -- there is no point in assigning a haplotype
  // to these, since they can be in any.
  //  vector<bool> unconstrained(numStrains, true); // this is too slow.
  memset(combined, '?', numStrains);

  // Computes number of defined alleles in block for each strain, and
  // initializes constrained strains in "combined".
  bool qMarks = qMarksInBlock(blockstart, blocksize); //check whether there are any unknowns for time saving later
  //int test;
  //char * combinedTest = (char*) malloc(numStrains);
  //memset(combinedTest, '?', numStrains);
  if (!qMarks)
  {
    return combinePatternNoQs(combined, blockstart, blocksize, haploLimit);
  }
  // numDefVec number of non-? entries in each column
  vector<int> numDefVec(numStrains);
  // Permutation of strain indices to sort columns by number of non-? entries.
  vector<int> strOrderVec(numStrains);

  // "pointers" to begining and end of block in snpVec, for "speed"
  vector<SNPInfo *>::iterator blockBeginIt = snpVec.begin() + blockstart;
  vector<SNPInfo *>::iterator blockEndIt = blockBeginIt + blocksize;
  //int total = blockEndIt - blockBeginIt;
  for (int str1 = 0; str1 < numStrains; str1++)
  {
    numDefVec[str1] = 0;
    // FIXME: Skip identical patterns here.
    // FIXME: I guess this could go by rows, which would be faster, too.
    for (vector<SNPInfo *>::iterator snpIt = blockBeginIt; snpIt != blockEndIt; snpIt++)
    {
      char *pat = (*snpIt)->pattern;
      char strchr = pat[str1];
      if (strchr != '?')
      {
        combined[str1] = '\000';
        numDefVec[str1]++;
      }
    }
  }

  // Sort a strOrderVec of indices 0..numStrains-1 so that they are in descending order
  // of numDefVec.  Leftmost strains are assigned first, later.
  vector<int>::iterator stoEnd = strOrderVec.end();
  int i = 0;
  for (vector<int>::iterator stoIt = strOrderVec.begin(); stoIt != stoEnd; stoIt++)
  {
    *stoIt = i++;
  }

  if (qMarks)
  {
    IndexComparator<int, greater<int>> idxCompare(&numDefVec);
    stable_sort(strOrderVec.begin(), strOrderVec.end(), idxCompare);
  }
  // Assign haplotypes to the strains.

  // This is a greedy method.  For each equivalence class index from 0 ... haploLimit
  // (or until done) choose leftmost strain with that equivalence class, assign same
  // equivalence class to all strains to the right of it.  Assign incompatible strains
  // to the right the next higher equivalence class (exit if haploLimit exceeded).
  // Invariants:  At end of loop doing EC i, all strains assigned to that class are
  // compatible nd will never change, and each is incompatible with the lower equivalence
  // classes, and eqclasse+1 is the maximum EC appearing in "combined" vector.
  // EC i+1 is incompatible with all lower equivalence classes.  In the inner
  // loop, strains that are processed satisfy this invariant.  The ones not done yet
  // satisfy the property for EC i-1.
  // In all the above, unconstrained strains have '?' in "combined" and are otherwise ignored.

  // FIXME: "done" check.

  char *merge = (char *)malloc(blocksize); //KSB what is the max? can we not malloc here? Why *in* for loop?
  vector<int>::iterator stoEnd1 = strOrderVec.end();

  int numHaplo = 0;
  for (int eqclass = 0; eqclass < haploLimit; eqclass++)
  {
    int firstStrainInEC = -1;
    memset(merge, '?', blocksize);
    // FIXME:  This doesn't need to start at beginning.  At least EQCLASS have already been
    // assigned.
    for (vector<int>::iterator stoIt1 = strOrderVec.begin() + eqclass; stoIt1 < stoEnd1; stoIt1++)
    { // not stoEnd, should be stoEnd1
      int str1 = *stoIt1;
      // skip strains that have already been assigned permanent ECs
      if (combined[str1] == eqclass)
      {
        // All other strains will be compared with this.
        if (firstStrainInEC < 0)
        {
          firstStrainInEC = str1;
          mergePattern(merge, blockstart, blocksize, str1);
        }
        else
        {
          // if incompatible, assign EC+1.  Otherwise, leave at current level.
          if (strainsAreCompatible(merge, blockstart, blocksize, str1))
          {
            // merge into combined column
            if (qMarks)
            {
              mergePattern(merge, blockstart, blocksize, str1);
            }
          }
          else
          {
            numHaplo = ++combined[str1] + 1; // WORRY: c++ ordering?
            if (numHaplo > haploLimit)
            {
              free(merge);

              /* if(!qMarks&&(test!=0&&!eqPatterns(combined, combinedTest))) {
                cout << "should be"; showPattern(combined);
                cout << endl;
                showPattern(combinedTest);
                cout << "from" << endl;
                showBlock(blockstart, blocksize);
                cout << "Returned " << test << " and limit is " << haploLimit << endl;
                cout << endl;
                    }
                    free(combinedTest);*/
              return 0;
            }
          }
        }
      }
    }
    // FIXME: This well-intended early exit did not take account of strvec reorder!
    // if (firstStrainInEC < 0 || firstStrainInEC >= numStrains-1) {
    if (firstStrainInEC < 0)
    {
      //free(merge);
      break; // Everything has been assigned.
    }
    //free(merge);
  }
  free(merge);

  // Makes colors ascend in order of strains.
  normalizePattern(combined);
  /*if(!qMarks&&(numHaplo!=test||!eqPatterns(combined, combinedTest))) {
    cout << "should be"; showPattern(combined);
    cout << endl;
    showPattern(combinedTest);
    cout << "from" << endl;
    cout << "Returned " << test << " instead of " << numHaplo << endl;
    cout << "limit is " << haploLimit << endl;
    cout << endl;
  }
  free(combinedTest);*/
  return numHaplo;
}

// Find highest-scoring block, starting at blockstart, only up to size maxSize
// and not exceeding haploLimit number of haplotypes.  Stores best block in by-ref first
// parameter.
void findBestBlock(HaploBlock &bestBlock, int blockstart, int haploLimit, int minSize, int maxSize = -1)
{
  if (traceFBB)
  {
    cout << "startBlock: " << blockstart << endl;
  }

  if (maxSize < 0)
  {
    maxSize = snpVec.size() - blockstart;
  }

  SNPInfo *pSNPInfo = snpVec[blockstart];

  // temporary block that we're scoring.
  HaploBlock tryBlock;
  tryBlock.start = blockstart;
  tryBlock.size = 1;
  tryBlock.score = .25; // 1 SNP/2^2 haplotypes
  memcpy(tryBlock.pattern, pSNPInfo->pattern, numStrains);

  bestBlock = tryBlock;

  int chrIdx = pSNPInfo->chrIdx;               // chromosome index of first SNP in block.
  int firstLoc = snpVec[blockstart]->position; // location on chromosome of first SNP.

  // look at sizes up to maxSize
  if (minSize < 2)
  {
    minSize = 2;
  }
  for (int blocksize = minSize; blocksize <= maxSize; blocksize++)
  {

    int numHaplo = combinePatterns(tryBlock.pattern, blockstart, blocksize, haploLimit);

    // break if too many haploTypes or chromosome changes.
    int curLoc = snpVec[blockstart + blocksize - 1]->position;
    if (numHaplo > 0 && chrIdx == snpVec[blockstart + blocksize - 1]->chrIdx && (curLoc - firstLoc <= tooLong))
    {
      tryBlock.size = blocksize;
      tryBlock.score = scoreBlock(blocksize, numHaplo);

      if (traceFBB)
      {
        cout << " Trying block size = " << blocksize << " -- score " << tryBlock.score << endl;
      }
      if (tryBlock.score >= bestBlock.score)
      {
        bestBlock = tryBlock; // tryBlock is the best so far.  If tied, take longer block.
      }
    }
    else
    {
      if (traceFBB)
      {
        cout << "Stop: numHaplo = " << numHaplo << ", span = " << curLoc - firstLoc << endl;
      }
      break;
    }
  }

  if (traceFBB)
  {
    cout << "Best block: " << bestBlock << endl;
  }
}

// Search for haplotype increase after power-of-two block
// This is functionally equivalent to findBestBlock if applied after a max score power-of-two block.
// minSize is the power-of-two block size -- we search after this.
// haploLimit in this case is the number of haplotypes in the power-of-two block.
void ffbBinSearch(HaploBlock &bestBlock, int blockstart, int haploLimit, int minSize, int maxSize = -1)
{
  // ssize is number of block sizes on each side of split that we
  // want to consider.
  int ssize = minSize / 2;

  int numHaplo = 0;

  char *combined = (char *)malloc(numStrains);

  // Look for haplotype increase before or after split.
  int split = minSize + ssize;

  int bestsize = minSize;

  SNPInfo *pSNPInfo = snpVec[blockstart];
  int chrIdx = pSNPInfo->chrIdx;
  int firstLoc = pSNPInfo->position;

  if (traceFBB)
  {
    cout << "  Binary search of power-of-2 block start = " << blockstart << ", size = " << minSize
         << " with " << haploLimit << " haplotypes." << endl;
    cout << "   Splitting at " << split << endl;
  }

  bool cached = false;
  // remember to check truncation. (blockend).
  while (true)
  {
    // is haplotype increase just before or after split?
    SNPInfo *pSplitSNPInfo = NULL;
    int splitLoc = -1;

    if (split <= maxSize)
    {
      // Make sure we're not going out of bounds before dereferencing
      pSplitSNPInfo = snpVec[blockstart + split - 1];
      splitLoc = pSplitSNPInfo->position;
      if (chrIdx != pSplitSNPInfo->chrIdx || (splitLoc - firstLoc > tooLong))
      {
        // Not allowed to use these SNPs.  Simulate haplotype increase.
        numHaplo = 0;
        if (traceFBB)
        {
          cout << "     split " << split << " is too far." << endl;
        }
      }
      else
      {
        numHaplo = combinePatterns(combined, blockstart, split, haploLimit);
        if (traceFBB)
        {
          cout << "   Combined patterns ";
          showPattern(combined);
          cout << ", numHaplo = " << numHaplo << endl;
        }
      }
    }
    else
    {
      // Not allowed to use these SNPs.  Simulate haplotype increase.
      numHaplo = 0;
      if (traceFBB)
      {
        cout << "     split " << split << " went past end of snpVec." << endl;
      }
    }

    if (numHaplo == 0)
    { // combinePatterns exceeded limit
      // best block index < split
      if (ssize == 1)
      {
        // split-1 is the best block
        bestsize = split - 1;
        if (traceFBB)
        {
          cout << "    best size is " << bestsize << endl;
        }
        break;
      }
      else
      {
        ssize /= 2;     // 1
        split -= ssize; // 5 -- want to see if 4 or 5 is best.
        if (traceFBB)
        {
          cout << "   Splitting at " << split << endl;
        }
      }
    }
    else
    {
      // best block index >= split.
      if (ssize == 1)
      {
        // split is the best block.
        bestsize = split;
        cached = true;
        if (traceFBB)
        {
          cout << "    best size is " << bestsize << endl;
        }
        break;
      }
      else
      {
        ssize /= 2;     // 1
        split += ssize; // 7 -- is best 6 or 7?
        if (traceFBB)
        {
          cout << "   Splitting at " << split << endl;
        }
      }
    }
  }

  double score = scoreBlock(bestsize, haploLimit);
  if (score >= bestBlock.score)
  {
    // Set out parameter bestBlock.
    bestBlock.start = blockstart;
    bestBlock.size = bestsize;
    bestBlock.score = score;
    if (cached)
    {
      memcpy(bestBlock.pattern, combined, numStrains);
    }
    else
    {
      combinePatterns(bestBlock.pattern, blockstart, bestsize, haploLimit);
    }
    if (traceFBB)
    {
      cout << "   New best Block " << bestBlock << endl;
    }
  }
  else
  {
    if (traceFBB)
    {
      cout << "  Locally but not globally best block." << endl;
    }
  }
  free(combined);
}

// Find highest-scoring block, starting at blockstart, only up to size maxSize
// and not exceeding haploLimit number of haplotypes.  Use size-doubling tricks
// to make it fast.
// WARNING: minSize is ignored.
void findBestBlockFast(HaploBlock &bestBlock, int blockstart, int haploLimit, int minSize, int maxSize = -1)
{
  if (traceFBB)
  {
    cout << "startBlock: " << blockstart << endl;
  }

  // Not sure why this was here, but it looks bogus.
  // compoundHaploBlocks.resize(snpVec.size());

  if (maxSize < 0)
  {
    maxSize = snpVec.size() - blockstart;
  }

  SNPInfo *pSNPInfo = snpVec[blockstart];

  // temporary block that we're scoring.
  HaploBlock tryBlock;
  tryBlock.start = blockstart;
  tryBlock.size = 1;
  tryBlock.score = .25; // 1 SNP, 2 haplotypes
  memcpy(tryBlock.pattern, pSNPInfo->pattern, numStrains);

  bestBlock = tryBlock;

  int chrIdx = pSNPInfo->chrIdx;               // chromosome index of first SNP in block.
  int firstLoc = snpVec[blockstart]->position; // location on chromosome of first SNP.

  // compute scores of power-of-two size blocks until there are too many haplotypes
  // or other stopping condition.  Also, compute the maximum score value seen.
  vector<HaploBlock> pow2blocks; // store info about the blocks for power-of-two sizes.
  double maxScore = 0.0;

  int prevNumHaplo = 0; // for monotonicity check

  for (int blocksize = 2; blocksize <= maxSize; blocksize *= 2)
  {
    int numHaplo = combinePatterns(tryBlock.pattern, blockstart, blocksize, haploLimit);

    // If coloring were exact, adding more SNPs to block would never reduce number of haplotypes.
    // However, coloring is approximate, and this algorithm can produce seriously suboptimal results
    // if the assumed monotonicity does not hold.  So this generates a warning.
    // TODO: Maybe there's a way to recover.
    if (false && numHaplo < prevNumHaplo && numHaplo > 0)
    {
      cout << "WARNING: Number of haplo blocks decreased with more SNPs!" << endl;
      cout << tryBlock << endl;
      cout << "Prev haplotypes: " << prevNumHaplo << ", cur haplotypes: " << numHaplo << endl;
    }
    prevNumHaplo = numHaplo;

    SNPInfo *pCurSNPInfo = snpVec[blockstart + blocksize - 1];
    int curLoc = pCurSNPInfo->position;
    // only proceed if this is an ok blocks (numHaplo == 0 means exceeded limit)
    if (numHaplo > 0 && chrIdx == pCurSNPInfo->chrIdx && (curLoc - firstLoc <= tooLong))
    {
      tryBlock.size = blocksize;
      tryBlock.score = scoreBlock(blocksize, numHaplo);

      if (tryBlock.score > maxScore)
      {
        maxScore = tryBlock.score;
      }

      if (traceFBB)
      {
        cout << " Trying block size = " << blocksize << " -- score " << tryBlock.score << endl;
      }

      pow2blocks.push_back(tryBlock);
    }
    else
    {
      if (traceFBB)
      {
        cout << "Stop: numHaplo = " << numHaplo << ", span = " << curLoc - firstLoc << endl;
      }
      break;
    }
  }

  // the true highest-scoring block has a size greater-equal to a maximum-scoring power-of-two-size
  // block and the less than the next power-of-two size.  If blocks are tied, search all of them.

  for (size_t i = 0; i < pow2blocks.size(); i++)
  {
    HaploBlock &pow2hb = pow2blocks[i];
    if (pow2hb.score == maxScore)
    {
      // It's a max scoring block.  Search for the truly best block after it and before
      // next one.
      // findBestBlock(bestBlock, pow2hb.start, haploLimit, pow2hb.size, maxSize);
      // derive haplo-limit from pow2hb (FIXME: could be faster?)
      int numHaplo = (int)(log2(pow2hb.size / pow2hb.score) + .5);
      ffbBinSearch(bestBlock, blockstart, numHaplo, pow2hb.size, maxSize);
    }
  }

  if (traceFBB)
  {
    cout << "Best block: " << bestBlock << endl;
  }
}

// tests whether block is too big: too many base pairs, off end of
// chromosome.  Precondition: Doesn't go past end of snpVec.
bool blockIsTooBig(int blockstart, int blocksize, int chrIdx, int firstLoc)
{
  SNPInfo *pCurSNPInfo = snpVec[blockstart + blocksize - 1];
  int curLoc = pCurSNPInfo->position;

  if (chrIdx != pCurSNPInfo->chrIdx)
  {
    if (traceFBB)
    {
      cout << "    End of chromosome" << endl;
    }
    return true;
  }
  else if (curLoc - firstLoc > tooLong)
  {
    if (traceFBB)
    {
      cout << "    Too many basepairs: " << curLoc - firstLoc << endl;
    }
    return true;
  }
  return false;
}

// Find the longest block beginning at blockstart of 1<= minsize <= size <= blocksize that has
// haploSize number of haplotypes.  blocksize should be legitimate -- it should
// not cross too large a gap between SNPs or cross over a chromosome boundary.
// Allocates and returns the haplolock with the right pattern.
// If it returns true, pHaploBlock is properly filled in with pattern, size, etc.
// Otherwise, pHaploBlock is not written.
bool fmbBinSearch(HaploBlock *pHaploBlock, int blockstart, int minSize, int blocksize, int haploSize, int chrIdx, int firstLoc)
{
  //cout << blockstart << "\t" << minSize << "\t" << blocksize << "\t" << haploSize << "\t" << chrIdx << "\t" << firstLoc << endl;
  if (blocksize < minSize)
  {
    return false;
  }

  // check if this size works.
  if (traceFBB)
  {
    cout << "fmbBinSearch:"
         << " blockstart: " << blockstart
         << ", minSize: " << minSize << ", blocksize: " << blocksize
         << ", haploSize: " << haploSize << endl;
  }

  // check end of chromosome, too many base pairs,
  if (blockIsTooBig(blockstart, blocksize, chrIdx, firstLoc))
  {
    return false;
  }

  char *tryPattern = (char *)malloc(numStrains);
  //memset(tryPattern, '?', numStrains); //combinePatterns does this
  int numHaplo = combinePatterns(tryPattern, blockstart, blocksize, haploSize);

  if (numHaplo == haploSize)
  {
    memcpy(pHaploBlock->pattern, tryPattern, numStrains);
    free(tryPattern);
    pHaploBlock->start = blockstart;
    pHaploBlock->size = blocksize;

    if (traceFBB)
    {
      cout << "   full block ok " << *pHaploBlock << endl;
    }
    return true;
  }

  free(tryPattern);

  if (minSize == blocksize)
  {
    // Block of haploSize does not exist.
    if (traceFBB)
    {
      cout << "   no block of size " << haploSize << "." << endl;
    }
    return false;
  }
  else if (minSize == blocksize - 1)
  {
    // only one possible size (minSize).
    return fmbBinSearch(pHaploBlock, blockstart, minSize, blocksize - 1, haploSize, chrIdx, firstLoc);
  }

  // Fall through: max size block <= blocksize-1, minSize <= blocksize-2

  int split = (blocksize + minSize) / 2; // minSize < split <= blocksize-1

  if (!fmbBinSearch(pHaploBlock, blockstart, minSize, split, haploSize, chrIdx, firstLoc))
  {
    // There is no block with haploSize (assuming monotonicity of combine).
    return false;
  }
  else if (pHaploBlock->size == split)
  {
    // full block up to split worked, so maybe we can find an even larger block.
    fmbBinSearch(pHaploBlock, blockstart, split + 1, blocksize - 1, haploSize, chrIdx, firstLoc);
  }

  // One of the two cases worked, so return true.  pHaploBlock was
  // set by one of them.
  return true;
}

// From a given start SNP find the longest block for a specified
// number of haplotypes (haploSize), while crossing too many
// base pairs or chromosome boundary.  Return a pointer to the HaploBlock,
// or NULL if no such block exists (when there is a jump of more than 1 in
// the number of haplotypes).
// Method:  Double block sizes until haploLimit is exceeded.  Then do binary
// search for longest block with specified number of haplotypes.  There may
// be no such block (number of haplotypes may increase by more than one when a
// SNP is added).
HaploBlock *findMaximalBlock(int blockstart, int haploSize, int minSize)
{
  if (traceFBB)
  {
    cout << "findMaximalBlock(" << blockstart << ", " << haploSize << ", " << minSize << ")" << endl;
  }

  int maxSize = snpVec.size() - blockstart;
  if (maxSize < minSize)
  {
    return NULL;
  }

  SNPInfo *pSNPInfo = snpVec[blockstart];

  int chrIdx = pSNPInfo->chrIdx;               // chromosome index of first SNP in block.
  int firstLoc = snpVec[blockstart]->position; // location on chromosome of first SNP.

  int prevNumHaplo = 0; // for monotonicity check
  int blocksize = 0;

  char *pattern = (char *)malloc(numStrains);
  memset(pattern, '?', numStrains);
  // Try progressively larger blocks until blocksize is definitely too big.
  for (blocksize = minSize; blocksize <= maxSize; blocksize *= 2)
  {
    if (traceFBB)
    {
      cout << "  Trying block size " << blocksize << endl;
    }

    if (blockIsTooBig(blockstart, blocksize, chrIdx, firstLoc))
    {
      break;
    }
    int numHaplo = combinePatterns(pattern, blockstart, blocksize, haploSize);
    // If coloring were exact, adding more SNPs to block would never reduce number of haplotypes.
    // However, coloring is approximate, and this algorithm can produce seriously suboptimal results
    // if the assumed monotonicity does not hold.  So this generates a warning.
    // TODO: Maybe there's a way to recover.
    // DISABLED:  Not easy to fix, doesn't seem to happen much.
    if (false && numHaplo < prevNumHaplo && numHaplo > 0)
    {
      cout << "WARNING: Number of haplo blocks decreased with more SNPs!" << endl;
      cout << pattern << endl;
      cout << "Prev haplotypes: " << prevNumHaplo << ", cur haplotypes: " << numHaplo << endl;
    }

    if (0 == numHaplo)
    {
      if (traceFBB)
      {
        cout << "    Too many haplotypes " << endl;
      }
      break;
    }
    prevNumHaplo = numHaplo;
  }
  free(pattern);

  // Loop only exits if blocksize is too big, so it's too big when we get here.
  // Make sure it won't go past end of SNPs.
  if (blocksize > maxSize)
  {
    blocksize = maxSize;
  }
  else if (blocksize <= minSize)
  {
    if (traceFBB)
    {
      cout << "   Cannot reduce block size due to minSize " << minSize << endl;
    }
    //return false;
    return NULL;
  }

  // Binary search for exact maximum block size.
  HaploBlock *pHaploBlock = new HaploBlock();
  pHaploBlock->start = 0;
  pHaploBlock->size = 0;
  int testBlock;
  if (blocksize == maxSize)
  {
    testBlock = blocksize;
  }
  else
  {
    testBlock = blocksize - 1;
  }

  if (!fmbBinSearch(pHaploBlock, blockstart, minSize, testBlock, haploSize, chrIdx, firstLoc))
  {
    // There is no block with the right number of haplotypes.
    delete pHaploBlock;
    pHaploBlock = NULL;
  }
  else
  {
    if (traceFBB)
    {
      cout << "Maximal block: " << *pHaploBlock << endl;
    }
  }
  //if(pHaploBlock != NULL){
  //showPattern(pHaploBlock->pattern);
  //cout << "\t" << pHaploBlock->start << "\t" << pHaploBlock->size  <<  endl;
  //}
  return pHaploBlock;
}

// Find best block starting at each SNP.
void findCompoundBlocks(int haploLimit)
{
  int blockstop = snpVec.size();
  compoundHaploBlocks.resize(snpVec.size());

  for (int blockstart = 0; blockstart < blockstop; blockstart++)
  {
    HaploBlock &hb = compoundHaploBlocks[blockstart];
    hb.start = blockstart;
    hb.size = 0;
    if (useFastFBB)
    {
      // optimization: If pattern is the same as previous SNP, don't
      // bother since block will end at the same place.  Furthermore,
      // if the previous block overlapped with a better preceding
      // block, the preceding block will have both SNPs, so it will
      // overlap this one, too.  The "tooLong" test is to deal with the case
      // where the pattern is the same as the previous SNP, but there is a big
      // gap that ends the previous block.  Also, look for chromosome transition.
      // FIXME: If this works, lift out of this IF
      if (0 == blockstart || snpVec[blockstart - 1]->pattern != snpVec[blockstart]->pattern || (snpVec[blockstart]->position - snpVec[blockstart - 1]->position) > tooLong || (snpVec[blockstart - 1]->chrIdx != snpVec[blockstart]->chrIdx))
      {
        findBestBlockFast(hb, blockstart, haploLimit, 0);
      }
      else if (traceFBB)
      {
        cout << "Skipping block because of identical patterns, blockstart = " << blockstart
             << ", pattern =";
        showPattern(snpVec[blockstart]->pattern);
        cout << endl;
      }
    }
    else
    {
      findBestBlock(hb, blockstart, haploLimit, 0);
    }
  }
}

// For each haplotype size, find the maximal length blocks.
// *** Add trace functions.
void findAllMaximalBlocks(int haploLimit, int minSNPBlocks)
{
  // Loop over haplotype numbers
  int minSize = 1; //minSNPBlocks;
  for (int haploSize = 2; haploSize <= haploLimit; haploSize++)
  {
    HaploBlock *prevSavedBlock = NULL;

    for (size_t blockstart = 0; blockstart < snpVec.size(); blockstart++)
    {
      // Only try sizes that aren't going to be subsumed!
      if (NULL != prevSavedBlock)
      {
        // Don't bother searching a block unless it will extend beyond previous saved block.
        minSize = prevSavedBlock->size - (blockstart - prevSavedBlock->start) + 1;
        if (minSize < 1)
        {
          minSize = 1;
        }
      }
      else
      {
        minSize = 1;
      }
      // minSize here is stride, that's how many step to move forward from prevSavedBlock
      HaploBlock *pHaploBlock = findMaximalBlock(blockstart, haploSize, minSize);

      if (NULL == pHaploBlock)
      {
      }
      else if (prevSavedBlock != NULL && pHaploBlock->size <= prevSavedBlock->size - (blockstart - prevSavedBlock->start))
      {
        // don't let the new block to be contained in the previous block,
        // previous end - current start is less than the size of the current block (size of the current block = current start - current end),
        // assumes previous start is always before current start - BY
        if (traceChooseBlocks)
        {
          cout << "Block " << *pHaploBlock << "  is subsumed by previous saved block " << *prevSavedBlock << endl;
        }
        delete pHaploBlock;
      }
      else
      {
        if (traceChooseBlocks)
        {
          cout << "Saving block " << *pHaploBlock << endl;
        }
        chosenHaploBlocks.push_back(pHaploBlock);
        prevSavedBlock = pHaploBlock;
      }
    }
  }

  // sort chosenHaploBlocks by chromosome and position.
  sort(chosenHaploBlocks.begin(), chosenHaploBlocks.end(), compareBlocksByPosition);
}

bool compareBySNP(HaploBlock *phb1, HaploBlock *phb2)
{
  SNPInfo *pSNP1 = snpVec[phb1->start];
  SNPInfo *pSNP2 = snpVec[phb2->start];
  return compareSNPs(pSNP1, pSNP2);
}

// Find best non-overlapping blocks.  Pick highest-score block, next highest-scoring
// that does not overlap, etc.
void chooseBlocks()
{
  // insert compoundHaploBlock pointers into priority queue.
  for (size_t blockstart = 0; blockstart < snpVec.size(); blockstart++)
  {
    sortedCompoundHaploBlocks.insert(&compoundHaploBlocks[blockstart]);
  }

  // Do all the blocks in decreasing score order.  For each candidate
  // block, scan SNPs in block to see if any is already marked as
  // "used", the current blocks overlaps a block that has already been
  // chosen.  In that case, rerun findBestBlock with start of current
  // block and maximumsize that goes up to but does not include the
  // used SNP.  Recompute the score and re-insert it with updated
  // score into the sortedCompoundHaploBlocks, to be tried again later
  // (since another block might have a higher score).  If no SNPs have
  // been used, mark them all as used and add this block to the chosen
  // blocks vector.

  chosenHaploBlocks.reserve(approxNumSNPs / 6);

  multiset<HaploBlock *, rCompareByScore>::iterator begit;
  while ((begit = sortedCompoundHaploBlocks.begin()) != sortedCompoundHaploBlocks.end())
  {

    begit = sortedCompoundHaploBlocks.begin();
    HaploBlock *phb = *begit;
    sortedCompoundHaploBlocks.erase(begit);

    if (traceChooseBlocks)
    {
      cout << "Candidate haplo block: " << *phb << endl;
    }

    if (phb->chosen)
    {
      cerr << "WARNING: chosen block on sortedCompoundHaploBlocks" << endl
           << "Block: " << *phb << endl;
    }

    // phb is a candidate for chosen block.
    int blocksize = phb->size;
    if (phb->size == 0)
    {
      // If size is 0, it was ruled out earlier, so just skip it.
      continue;
    }
    int blockstart = phb->start;
    bool choose = true; // choose unless it overlaps another chosen block.

    // check whether one of phb's blocks has already been chosen, and truncate if necessary.
    for (int j = 0; j < blocksize; j++)
    {

      SNPInfo *pSNPInfo = snpVec[blockstart + j]; // overlapped block.

      if (pSNPInfo->used)
      {
        if (j > 0)
        {
          // phb overlaps an already chosen block.  Truncate and toss it back.
          double oldscore = phb->score;
          // j is max size of block.  If start=10, j=1, SNP 11 is used, then
          // our max size block is j=1.
          if (traceChooseBlocks)
          {
            cout << "Truncating block to size " << j << ":" << endl
                 << *phb << endl;
          }
          if (useFastFBB)
          {
            findBestBlockFast(compoundHaploBlocks[blockstart], blockstart, haploLimit, 0, j);
          }
          else
          {
            findBestBlock(compoundHaploBlocks[blockstart], blockstart, haploLimit, 0, j);
          }
          // Score of block has changed, so reinsert it
          // Invariant: no chosen blocks in sortedCompoundHaploBlocks.
          if (phb->chosen)
          {
            cout << "Reinserting chosen phb???" << endl;
            cout << *phb << endl;
          }
          if (traceChooseBlocks)
          {
            cout << "Reinserting block: " << *phb << endl;
          }
          sortedCompoundHaploBlocks.insert(phb);
          if (phb->score > oldscore)
          {
            cerr << "WARNING: Truncated block has higher score than original score of "
                 << oldscore << endl
                 << *phb << endl;
          }
        }
        else
        {
          // our first block is overlapped by a preceding block.  Just give up.
          if (traceChooseBlocks)
          {
            cout << "  ... overlapped " << endl;
          }
        }
        choose = false;
        break;
      }
    }

    // None of phb's blocks were already chosen, so choose phb.
    if (choose)
    {
      // phb did not overlap a chosen block, so phb should be chosen.
      chosenHaploBlocks.push_back(phb);
      phb->chosen = true;
      // mark overlapped blocks.
      for (int j = 0; j < blocksize; j++)
      {
        SNPInfo *pSNPInfo = snpVec[blockstart + j];
        if (pSNPInfo->used)
        {
          cerr << "WARNING: About to mark use SNPInfo as used." << endl
               << "SNPInfo: " << *pSNPInfo << endl;
        }
        pSNPInfo->used = true;
        // cout << "Marked overlapped: " << *pSNPInfo << endl;
      }
      if (traceChooseBlocks)
      {
        cout << "Chosen: " << *phb << endl;
      }
    }
  }
  // sort chosen blocks by chromosome/position.
  sort(chosenHaploBlocks.begin(), chosenHaploBlocks.end(), compareBySNP);
}

// print a block and the SNPs in it.
void showBlockSNPs(HaploBlock &hb)
{
  cout << hb;
  for (int i = hb.start; i < hb.start + hb.size; i++)
  {
    cout << "  " << i << ":\t" << *snpVec[i] << endl;
  }
}

void showBlocks(vector<HaploBlock> &hbv)
{
  for (size_t i = 0; i < hbv.size(); i++)
  {
    showBlockSNPs(hbv[i]);
  }
}

// Show a vector of pointers to blocks, instead of blocks. (Badly named.)
void showBlocksSNPs(vector<HaploBlock *> &hbs)
{
  for (size_t i = 0; i < hbs.size(); i++)
  {
    showBlockSNPs(*hbs[i]);
  }
}

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

  friend ostream &operator<<(ostream &os, const StrainDisp &sd);
};

ostream &operator<<(ostream &os, const StrainDisp &sd)
{
  os << "StrainDisp " << strainAbbrevs.eltOf(sd.strIdx) << "(" << sd.haplotype << ")"
     << ": strIdx = " << sd.strIdx
     << ", numStrainsInSameHaplotype = " << sd.numStrainsInSameHaplotype
     << ", minStrainInSameHaplotype = " << sd.minStrainInSameHaplotype
     << ", numSNPsDefined = " << sd.numSNPsDefined
     << flush;
  return os;
};

// comparison function to get in right order.
// haplotypes in decreasing order of number of strains.
// w/in haplotype, put most defined ones first?
//   haplotype -> straincount array.
//   strIdx (order in strain_index file) to break ties.
//  FIXME: Why can't I used references for sd1, sd2?
bool compareStrainDisps(StrainDisp sd1, StrainDisp sd2)
{
  if (sd1.numStrainsInSameHaplotype != sd2.numStrainsInSameHaplotype)
  {
    return (sd1.numStrainsInSameHaplotype > sd2.numStrainsInSameHaplotype);
  }
  else if (sd1.minStrainInSameHaplotype != sd2.minStrainInSameHaplotype)
  {
    // I'd like to order haplotypes by minimum strIdx,  but I'm tired.
    return (sd1.minStrainInSameHaplotype < sd2.minStrainInSameHaplotype);
  }
  // STRANGE: sort would trash vector elements when this incorrectly said ">" instead
  // of ">".  Hard to imagine why (might not have been algebraically a comparison?)
  else if (sd1.numSNPsDefined != sd2.numSNPsDefined)
  {
    return (sd1.numSNPsDefined > sd2.numSNPsDefined);
  }
  else
  {
    return (sd1.strIdx <= sd2.strIdx);
  }
}

// Write out html for one block.
void writeBlock(ostream &hs, int blkIdx, int haploLimit)
{
  HaploBlock *hb = chosenHaploBlocks[blkIdx];
  char *pattern = hb->pattern;
  int blockstart = hb->start;
  int blocksize = hb->size;

  int numHaplo = countHaplotypes(pattern);

  // vector mapping haplotype to the number of strains in that haplotype
  vector<int> numHaplotypeStrains(numHaplo, 0);
  // minimum strain index in each haplotype.
  vector<int> minHaplotypeStrain(numHaplo, numStrains);

  // vector of strain displays
  vector<StrainDisp> strainDispVec(numStrains);

  // Set up StrainDisp objects in strainDispVec
  for (int strIdx = 0; strIdx < numStrains; strIdx++)
  {
    int haplotype = pattern[strIdx];
    if (haplotype != '?')
    {
      numHaplotypeStrains[haplotype]++;
      if (minHaplotypeStrain[haplotype] > strIdx)
      {
        minHaplotypeStrain[haplotype] = strIdx;
      }
    }
    // strainDispVec.push_back(StrainDisp());
    // StrainDisp & strainDisp = strainDispVec.back();
    StrainDisp &strainDisp = strainDispVec[strIdx];
    strainDisp.strIdx = strIdx;
    strainDisp.haplotype = haplotype;

    vector<SNPInfo *>::iterator blockBeginIt = snpVec.begin() + blockstart;
    vector<SNPInfo *>::iterator blockEndIt = blockBeginIt + blocksize;
    for (vector<SNPInfo *>::iterator snpIt = blockBeginIt; snpIt != blockEndIt; snpIt++)
    {
      SNPInfo *pSNPInfo = *snpIt;
      if (pSNPInfo->alleles[strIdx] != '?')
      {
        strainDisp.numSNPsDefined++;
      }
    }
  }

  // Now that we know number of strains in each haplotype, set those
  // in StrainDisps
  for (size_t strIdx = 0; strIdx < strainDispVec.size(); strIdx++)
  {
    StrainDisp &strainDisp = strainDispVec[strIdx];
    int haplotype = pattern[strIdx];
    if (haplotype != '?')
    {
      strainDisp.numStrainsInSameHaplotype = numHaplotypeStrains[haplotype];
      strainDisp.minStrainInSameHaplotype = minHaplotypeStrain[haplotype];
    }
  }

  //  for (int j = 0; j < strainDispVec.size(); j++) {
  // cout << "StrainDisp " << j << " strIdx = " << strainDispVec[j].strIdx << endl;
  // }

  sort(strainDispVec.begin(), strainDispVec.end(), compareStrainDisps);

  if (traceWriteHTML)
  {
    cout << "Block " << *hb << endl;
    for (size_t j = 0; j < strainDispVec.size(); j++)
    {
      cout << "  " << strainDispVec[j] << endl;
    }
  }

  // begin table
  hs << "<HR><TABLE BORDER=0 CELLSPACING=3>" << endl;
  hs << "<TR><TH><FONT SIZE=1><A name=block_" << blkIdx << ">block "
     << blkIdx << "&nbsp;&nbsp;&nbsp;&nbsp;</A></FONT></TH>";
  hs << "<TH><FONT COLOR=white>XXX</FONT></TH><TH><FONT COLOR=white>XXXX</FONT></TH>";

  int prevHaplotype = 0;

  // print strain names row.
  for (size_t i = 0; i < strainDispVec.size(); i++)
  {
    StrainDisp &strainDisp = strainDispVec[i];
    int strIdx = strainDisp.strIdx;
    int haplotype = pattern[strIdx];

    if (haplotype != prevHaplotype)
    {
      // extra space between haplotypes.
      hs << "<TD><FONT COLOR=white SIZE=1>X</FONT></TD>";
      prevHaplotype = haplotype;
    }

    // print abbreviated strain name.  Gray if there are any undefined SNPs for this strain.
    hs << "<TH>";
    if (strainDisp.numSNPsDefined < blocksize)
    {
      hs << "<FONT COLOR=gray>";
    }
    hs << strainAbbrevs.eltOf(strainDisp.strIdx);
    if (strainDisp.numSNPsDefined < blocksize)
    {
      hs << "</FONT>";
    }
    hs << "</TH>";
  }
  hs << "</TR>" << endl;

  // print the rows for SNPs
  for (int snpIdx = hb->start; snpIdx < hb->start + hb->size; snpIdx++)
  {
    SNPInfo *pSNPInfo = snpVec[snpIdx];
    char *alleles = pSNPInfo->alleles;

    hs << "<TR><TD><FONT COLOR=red SIZE=1></FONT></TD><TD><FONT COLOR=green SIZE=1>"
       << snpIdx << "</FONT></TD><TD></TD>";

    prevHaplotype = 0;

    char majorAllele, minorAllele;
    // major allele is blue, minor allele is yellow.
    findMajorMinorAllele(alleles, majorAllele, minorAllele);

    for (size_t i = 0; i < strainDispVec.size(); i++)
    {
      StrainDisp &strainDisp = strainDispVec[i];
      int strIdx = strainDisp.strIdx;
      int haplotype = pattern[strIdx];

      if (haplotype != prevHaplotype)
      {
        // extra space between haplotypes.
        hs << "<TD><FONT COLOR=white SIZE=1>X</FONT></TD>";
        prevHaplotype = haplotype;
      }

      char allele = alleles[strIdx];
      hs << "<TD BGCOLOR=";

      if ('?' == allele)
      {
        hs << "white><FONT COLOR=white SIZE=1>0";
      }
      else if (allele == majorAllele)
      {
        hs << "blue><FONT COLOR=blue SIZE=1>2";
      }
      else if (allele == minorAllele)
      {
        hs << "yellow><FONT COLOR=yellow SIZE=1>1";
      }
      hs << "</FONT></TD>";
    }
    // very end of line for SNP
    hs << "<TD><FONT SIZE=1>" << pSNPInfo->position << "</FONT></TD>" << endl;
  }

  // end table
  hs << "</TABLE>" << endl;
}

// Writes HTML pages for a chromosome.  Since SNPs for multiple
// chromosomes may be in same vector and chosenBlocks may have the
// same property, this takes blkIdx, which is the index of the
// chosenBlock where the chromosome starts, and it returns the
// beginning of the next chromosome (or size+1 if at end). Returns
// first block of next chromosome (chosenHaploBlocks.size() when
// done).

int writeChromosome(const char *dirname, size_t blkIdx)
{
  int chrIdx = snpVec[chosenHaploBlocks[blkIdx]->start]->chrIdx;

  string chrName = chromosomes.eltOf(chrIdx);

  string fname(dirname + chrName + ".html");
  // string fname("/hd/dill/users/dill/WWW/haplo/" + chrName + ".html");
  // string fname("haplo/" + chrName + ".html");

  ofstream hs(fname.c_str());
  if (!hs.is_open())
  {
    cerr << "Open of file \"" << fname << "\" failed: ";
    perror("");
    exit(1);
  }

  hs << "<HTML>" << endl
     << "<HEAD><TITLE>Chromosome " << chrName << " Haplotype</TITLE></HEAD>" << endl
     << "<BODY BGCOLOR=white VLINK=red TEXT=black LINK=red ALINK=green>" << endl
     << "<CENTER>" << endl
     << "<FONT SIZE=8 FACE=Arial>Chromosome " << chrName << " </FONT>" << endl
     << "<P><HR><P>" << endl;

  // Using a parameter as a loop variable.  Is that bad?
  for (;
       blkIdx < chosenHaploBlocks.size() && (chrIdx == snpVec[chosenHaploBlocks[blkIdx]->start]->chrIdx);
       blkIdx++)
  {
    writeBlock(hs, blkIdx, haploLimit);
  }

  hs << "<HR></CENTER></BODY></HTML>" << endl;

  hs.close();

  return blkIdx;
}

void writeHTML(char *dirName)
{
  for (size_t blkIdx = 0; blkIdx < chosenHaploBlocks.size(); blkIdx = writeChromosome(dirName, blkIdx))
    ;
}

// write out list of gene names for a single block.
void writeBlockGeneNames(ofstream &os, HaploBlock *pHB)
{
  map<string, string> geneIsCodingMap;
  // Look through all SNPs in block to find overlapping genes, record whether they
  // include non-synonymous changes.
  size_t snpEnd = pHB->start + pHB->size;
  for (size_t snpIdx = pHB->start; snpIdx < snpEnd; snpIdx++)
  {
    SNPInfo *pSNPInfo = snpVec[snpIdx];
    map<string, string> &geneCodonMap = pSNPInfo->geneCodonMap;
    //ofstream debug_log;
    //debug_log.open("debug.log",ios::app);

    for (map<string, string>::iterator gcit = geneCodonMap.begin(); gcit != geneCodonMap.end(); gcit++)
    {
      //cout << (*gcit).first << "\t" << gcit->first << "\t" << (*gcit).second << "\t" << pSNPInfo->position << endl;
      if (geneIsCodingMap[(*gcit).first] == "")
      {
        //debug_log << "The first time: " << gcit->first << " " << gcit->second << endl;
        if ((*gcit).second.find("NON_SYNONYMOUS_CODING") != string::npos)
        {
          geneIsCodingMap[(*gcit).first] = "1";
        }
        else if ((*gcit).second.find("<") != string::npos || (*gcit).second.find("SPLICE_SITE") != string::npos)
        {
          geneIsCodingMap[(*gcit).first] = (*gcit).second;
        }
        else
        {
          geneIsCodingMap[(*gcit).first] = "0";
        }
      }
      else
      {
        //debug_log << "NOT the first time: " << gcit->first << " " << gcit->second << endl;
        if ((*gcit).second.find("NON_SYNONYMOUS_CODING") != string::npos)
        {
          if (geneIsCodingMap[(*gcit).first] == "0")
            geneIsCodingMap[(*gcit).first] = "1";
          //continue;
        }
        if ((*gcit).second.find("<") != string::npos || (*gcit).second.find("SPLICE_SITE") != string::npos)
        {
          //debug_log << "splice site in" << gcit->first << endl;
          if (geneIsCodingMap[(*gcit).first] == "0" || geneIsCodingMap[(*gcit).first] == "1")
          {
            geneIsCodingMap[(*gcit).first] = (*gcit).second;
            //debug_log << "changing from 0/1: gene name: " << gcit->first << "descriptions: " << geneIsCodingMap[gcit->first] << endl;
          }
          else
          {
            geneIsCodingMap[(*gcit).first] = geneIsCodingMap[(*gcit).first] + "!" + (*gcit).second;
            //debug_log << "gene name: " << gcit->first << "descriptions: " << geneIsCodingMap[gcit->first] << endl;
          }
        }
      }
    }
    //debug_log.close();
  }
  for (map<string, string>::iterator gicit = geneIsCodingMap.begin(); gicit != geneIsCodingMap.end(); gicit++)
  {
    os << "\t" << (*gicit).first << "\t" << (*gicit).second;
  }
}

// Write block summary for a single chromosome.
int writeChrBlockSummary(ofstream &os, size_t blkIdx, int minBlockSNPs)
{
  int chrIdx = snpVec[chosenHaploBlocks[blkIdx]->start]->chrIdx;

  string chrName = chromosomes.eltOf(chrIdx);

  // Using a parameter as a loop variable.  Is that bad?
  for (;
       blkIdx < chosenHaploBlocks.size() && (chrIdx == snpVec[chosenHaploBlocks[blkIdx]->start]->chrIdx);
       blkIdx++)
  {

    HaploBlock *pHB = chosenHaploBlocks[blkIdx];
    // FIXME: Write out all blocks, regardless of size. Filter out small blocks in phmap, not here.

    if (pHB->size >= minBlockSNPs)
    {

      SNPInfo *pFirstSNP = snpVec[pHB->start];
      SNPInfo *pLastSNP = snpVec[pHB->start + pHB->size - 1];

      int qs = 0;
      for (int i = pHB->start; i < (pHB->start + pHB->size); i++)
      {
        if (snpVec[i]->qMarks)
        {
          qs = 1;
          break;
        }
      }
      os << chrName << "\t" << blkIdx << "\t" << pHB->start << "\t"
         << pHB->size << "\t" << pFirstSNP->position << "\t"
         << pLastSNP->position << "\t";
      showPattern(os, pHB->pattern);
      writeBlockGeneNames(os, pHB);
      //os << "\t" << qs;
      os << endl;
    }
  }

  return blkIdx;
}

int acquireLock()
{
  // int result;
  //*** Disable locking
  return 0;
  // while (result = open("output.lck", O_WRONLY | O_CREAT | O_EXCL, S_IRUSR|S_IWUSR) < 0 && errno==EEXIST) {
  //   sleep(1);
  // }
  // if (result<0) { //oh noes
  //   cerr << "error in acquiring file lock to write output" << endl;
  //   perror("");
  //   exit(1);
  // }
  // return result;
}

void releaseLock(int fid)
{
  // close(fid);
  return;
  //  unlink("output.lck");
}

void writeBlockSummary(char *fileName, int minBlockSNPs)
{
  int fid = acquireLock();
  ofstream os(fileName);
  if (!os.is_open())
  {
    cerr << "Open of file \"" << fileName << "\" failed: ";
    perror("");
    exit(1);
  }

  // Header made concatenation of individual chromosome results more difficult.
  // os << "Chr\tBlock\tStart\tSize\tChrBeg\tChrEnd\tPattern" << endl;
  // FIXME: Since nhaploblocks is only applied to one chromosome at a time, writing out the
  // chromosomes separately as this does is probably not useful.
  for (size_t blkIdx = 0; blkIdx < chosenHaploBlocks.size(); blkIdx = writeChrBlockSummary(os, blkIdx, minBlockSNPs))
    ;
  releaseLock(fid);
}

// Writes SNPs that were actually used in the block, in order of chromosome/position.
// This is specific to the set of strains and minBlockSNPs, because the good SNPs depend on the strains.
// This may include some SNPs in blocks that were discarded for being too small.
// FORMAT:  chrom	Pos	SNPID	Pattern		geneName1	codon1	geneName2 ...
void writeBlockSNPs(char *fname, int minBlockSNPs)
{
  int fid = acquireLock();
  ofstream snpOut(fname);
  if (!snpOut.is_open())
  {
    cerr << "Open of file \"" << fname << "\" failed: ";
    perror("");
    exit(1);
  }

  for (size_t snpIdx = 0; snpIdx < snpVec.size(); snpIdx++)
  {
    SNPInfo *pSNPInfo = snpVec[snpIdx];
    string chrName = chromosomes.eltOf(pSNPInfo->chrIdx);
    snpOut << chrName << "\t" << pSNPInfo->position << "\t"
           << pSNPInfo->name << "\t";
    showPattern(snpOut, pSNPInfo->pattern);
    map<string, string> &geneMap = pSNPInfo->geneCodonMap;
    for (map<string, string>::iterator gcit = geneMap.begin(); gcit != geneMap.end(); gcit++)
    {
      snpOut << "\t" << (*gcit).first << "\t" << (*gcit).second;
    }
    snpOut << endl;
  }
  releaseLock(fid);
}

bool isCompatiblePattern(char *p1, char *p2)
{
  bool result = true;
  for (int str = 0; str < numStrains; str++)
  {
    if (p1[str] != p2[str] && p1[str] != '?' && p2[str] != '?')
    {
      return (false);
    }
  }
  return result;
}

// update stats after we have a chosenHaploBlocks vector (in case we do it after
// every chromosome).
void updateStats()
{
  totalBlocks += chosenHaploBlocks.size();
  totalSNPs += snpVec.size();
  for (size_t blkIdx = 0; blkIdx < chosenHaploBlocks.size(); blkIdx++)
  {
    HaploBlock *pHB = chosenHaploBlocks[blkIdx];
    totalHaploTypes += countHaplotypes(pHB->pattern);
    if (pHB->size <= 3)
    {
      // small block;
      numSmallBlocks++;
      totalSNPsInSmallBlocks += pHB->size;
    }
  }
}

void reportStats()
{
  cout << "total blocks:\t" << totalBlocks << endl;
  cout << "total SNPs:\t" << totalSNPs << endl;
  cout << "avg haplotypes/block:\t" << ((double)totalHaploTypes) / ((double)totalBlocks)
       << endl;
  cout << "number of small blocks (1, 2, or 3 SNPs):\t" << numSmallBlocks << endl;
  cout << "total SNPs in small blocks:\t" << totalSNPsInSmallBlocks << endl;
}

class Options
{

public:
  bool non_overlapping;
  char *strainsFileName;
  char *allelesFileName;
  char *genesFileName;
  char *blocksFileName;
  char *blockSNPsFileName;
  int minBlockSNPs;

  // constructor
  Options() : strainsFileName(NULL), blocksFileName(NULL){};
};

Options *parseOptions(int argc, char **argv)
{
  int c;

  Options *opts = new Options();

  static struct option long_options[] = {
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
      {"version", no_argument, 0, 0},
      {0, 0, 0, 0}};

  // default values
  opts->minBlockSNPs = 4;

  const char *usage = "usage:\n eblocks\n"
                      "    -v --verbose\n"
                      "    -n --non_overlapping\n"
                      "    -s --strains_file <name of strains file>\n"
                      "    -a --alleles_file <name of alleles file>\n"
                      "    -g --genes_file <name of SNPs-to-genes file>\n"
                      "    -o --blocks_file  <output file with haplotype blocks>\n"
                      "    -p --blockSNPs_file  <output file with SNPs for blocks>\n"
                      "    -m --min_block_SNPs\n"
                      "    -h --help\n"
                      "    --version\n";

  while (1)
  {

    int option_index = 0;
    c = getopt_long(argc, argv, "hvna:g:s:o:p:m:", long_options, &option_index);

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
    case 0: /* long option without a short arg */
    {
      if (strcmp("version", long_options[option_index].name) == 0)
      {
        std::cout << "eblock version: " << __VERSION__ << std::endl;
      }
      exit(1);
    }
    default:
      abort();
    }
  }

  // isn't there some way to do this automatically based on above info?
  if (NULL == opts->strainsFileName)
  {
    cerr << "Required arg missing: strains file name (-s)" << endl;
    cerr << usage << endl;
    exit(1);
  }
  else if (NULL == opts->allelesFileName)
  {
    cerr << "Required arg missing: alleles file name (-a)" << endl;
    cerr << usage << endl;
    exit(1);
  }
  else if (NULL == opts->blocksFileName)
  {
    cerr << "Required arg missing: blocks file name (-o)" << endl;
    cerr << usage << endl;
    exit(1);
  }
  else if (NULL == opts->blockSNPsFileName)
  {
    cerr << "Required arg missing: block SNPs file name (-p)" << endl;
    cerr << usage << endl;
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

// Main for nhaploblocks
int main(int argc, char **argv)
{
  Options *opts = parseOptions(argc, argv);

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
  string chr = chromosomes.eltOf(0);
  // if no gene names file specified, don't read it, omit gene names.
  if (opts->genesFileName)
  {
    beginPhase("reading SNP gene names file");
    // FIXME: make this an option.
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
  beginPhase("writing block SNPs");
  writeBlockSNPs(opts->blockSNPsFileName, opts->minBlockSNPs);
  endPhase();

  if (verbose)
  {
    updateStats();
    reportStats();
  }

  return 0;
}

// This just compacts the perlegen tables.
// KEEP THIS -- needed to compact perlegen tables
int main3(int argc, char **argv)
{

  //  cout << "Just doing chromosome 1, for testing" << endl;
  //  beginPhase("reading chr_index.txt");
  //  readDynumList("PERLEGEN/chr_index.txt", chromosomes);
  //  endPhase();

  beginPhase("reading strains");
  // read all perlegen strains.
  readStrains((char *)"all_strains.txt");
  endPhase();

  numStrains = strainAbbrevs.size();

  // perlegen_b36_snp_vs_mmgene_091208.unl
  beginPhase("reading perlegen_b36_snp_vs_mmgene_091208.unl");
  readPerlegenSNPvsmgene((char *)"perlegen_b36_snp_vs_mmgene_091208.unl");
  // readPerlegenSNPvsmgene("perlegen_test_chr.unl");
  endPhase();

  beginPhase("reading perlegen_b36_strain_091208.unl");
  readPerlegenAlleleInfo((char *)"perlegen_b36_strain_091208.unl");
  //readPerlegenAlleleInfo("perlegen_test_snps.unl");
  endPhase();

  beginPhase("writing perlegen_b36_compact.txt");
  writeAlleleInfoCompact((char *)"perlegen_b36_compact.txt");
  // readPerlegenAlleleInfo("perlegen_test_snps.unl");
  endPhase();

  return 0;
}

// compact the roche tables
// KEEP THIS -- you need it to build the compact tables.
int main4(int argc, char **argv)
{

  beginPhase("reading chr_index.txt");
  // readDynumList("/home/dill/bin/data/chr_index.txt", chromosomes);
  // readDynumList("/gourd/dill/bin/data/chr_index.txt", chromosomes);
  readDynumList((char *)"chr_index.txt", chromosomes);
  endPhase();
  // this has the strains we actually want to use.
  // readDynumList("strain_index.txt", relevantStrains);

  //  beginPhase("reading h2_strain_index.txt");
  // readStrains("h2_strain_index.txt");
  beginPhase("reading strains file");
  readStrains((char *)"roche_all_strains.txt");
  endPhase();

  numStrains = strainAbbrevs.size();
  minDefined = (numStrains + 1) / 2; // half of strains must be defined (rounded up).

  beginPhase("reading chr_info_perl.txt");
  // readChromosomeInfo("/gourd/dill/bin/data/chr_info_perl.txt");
  readChromosomeInfo((char *)"chr_info_perl.txt");
  endPhase();

  beginPhase("reading allele_info_perl.txt");
  // readAlleleInfo("/gourd/dill/bin/data/allele_info_perl.txt");
  readAlleleInfo((char *)"allele_info_perl.txt");
  // readAlleleInfo("allele_info_test.txt");
  endPhase();

  beginPhase("writing roche_compact.txt");
  writeAlleleInfoCompact((char *)"roche_compact.txt");
  // readPerlegenAlleleInfo("perlegen_test_snps.unl");
  endPhase();

  return 0;
}
