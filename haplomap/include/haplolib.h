// -*-c++-*-
#ifndef HAPLOLIB_H
#define HAPLOLIB_H
// Various useful functions from nhaploblocks
//#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale> // for numbers with commas in them.
#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cmath>
#include <time.h>
#include <unordered_set>
//#include "utils.h"
#include "dynum.h"
#include "indexcomparator.h"
#include "ColumnReader.h"

//using namespace std;

// Copy a pattern string.
char *dupPattern(char *pattern);

std::ostream &showPattern(const char *pattern);
std::ostream &showPattern(const char *pattern, size_t len);
std::ostream &showPattern(std::ostream &os, const char *pattern);
std::ostream &showPattern(std::ostream &os, const char *pattern, size_t len);
void readStrains(char *fname);

// Progress messages.
void beginPhase();
void beginPhase(const char *msg);
void endPhase();
void endPhase(const char *msg, std::string chr);


// Everything we need to know about a SNP.
class SNPInfo
{
    char *initAlleles();
    // has to be called after initAlleles.
    char *initPattern();

public:
    friend std::ostream &operator<<(std::ostream &os, const SNPInfo &s);

    // default constructor
    SNPInfo();
    // constructor for stuff in chr_info_perl file.
    //  SNPInfo(string n, int ci, int pos) : name(n), chrIdx(ci), position(pos) { alleles = initAlleles(); };
    SNPInfo(std::string n, int ci, int pos);

    // copy constructor
    SNPInfo(SNPInfo const &snpInfo);
    // destructor
    ~SNPInfo();
    // Overload assignment to memcpy, so we don't have to worry about string copying.
    SNPInfo &operator=(const SNPInfo &si);

    std::string name;
    int chrIdx;   // chromosome index
    int position; // location on chromosome
    char ref; // reference allele char
    // string of alleles "ACGTD?" "D" is deletion, "?" is unspecified.
    char *alleles;
    // string of haplotype numbers (not printable).
    char *pattern;
    std::map<std::string, std::string> geneCodonMap; // Map gene names to codons/annotations
    bool frozen;                      // when true, pattern is uniquified and should not be changed.
    bool used;                        // Is this SNP in a chosen block?
    bool qMarks;                      //Are there any unkowns in this SNP?

    char getAllele(int i);
    void setAllele(int i, char a);
    void setReference(char c);
    char getReference(void);

    // Convert allele string (e.g., "?ACCAT?") to a pattern (e.g., "?01102?").
    // This is an updated version of allelesToPattern. 
    // only used for writting output pattern, explicity set (ref: 0, alternate: 1-3, unknown: ?)
    // old version ? we don't know 0 is ref or alternate.
    std::string allelesToPattern(void);

};

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

    void init();

    // default constructor
    HaploBlock();

    // copy constructor
    HaploBlock(HaploBlock const &hb);

    // Overload assignment to memcpy, so we don't have to worry about string copying.
    HaploBlock &operator=(const HaploBlock &hb);

    // destructor
    ~HaploBlock();

    friend std::ostream &operator<<(std::ostream &os, const HaploBlock &s);
};



struct rCompareByScore
{
    bool operator()(const HaploBlock *phb1, const HaploBlock *phb2) const
    {
        return phb1->score > phb2->score;
    }
};


// Test patterns for exact equality
bool eqPatterns(const char *pat1, const char *pat2);

// hash function -- standard char* hash fails because of null termination
struct PatternHash
{
    size_t operator()(const char *s) const;
};

// Hash table for pattern uniquification
struct PatternEq
{
    bool operator()(const char *pat1, const char *pat2) const
    {
        return eqPatterns(pat1, pat2);
    }
};


// for debugging
inline void showBlock(const HaploBlock &hb)
{
    std::cout << hb << std::endl;
}

void showBlocks(std::vector<HaploBlock *> &hbv);


template <typename T>
decltype(std::bind(&T::value_type::second, std::placeholders::_1)) select2nd() {
    return std::bind(&T::value_type::second, std::placeholders::_1);
}
// with parameters
//template <typename T>
//decltype(std::bind(&T::value_type::second, std::placeholders::_1)) select2nd(T m) {
//    return std::bind(&T::value_type::second, std::placeholders::_1);
//}

//typedef hash_set<char*, PatternHash, PatternEq> PatternSet;
typedef std::unordered_set<char *, PatternHash, PatternEq> PatternSet;
// FIXME: tasteless global variables (maybe verbose is ok).

extern int numStrains;
extern bool verbose; // in haploLib

extern Dynum<std::string> relevantStrains;
extern Dynum<std::string> strainAbbrevs;
extern PatternSet patternUniqueTable;
// Map chromosome names to numerical indices.
extern Dynum<std::string> chromosomes;
// vector of all good SNPs in order of chromosome, position.
extern std::vector<SNPInfo *> snpVec;

#endif