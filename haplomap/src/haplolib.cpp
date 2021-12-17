// -*-c++-*-

// Various useful functions from nhaploblocks

#include "haplolib.h"
#include <unistd.h>

static time_t start_time; // For time printing.  Should be in an object.

int numStrains = -1;
Dynum<std::string> relevantStrains;
Dynum<std::string> strainAbbrevs;

// Map chromosome names to numerical indices.
Dynum<std::string> chromosomes;
// vector of all good SNPs in order of chromosome, position.
std::vector<SNPInfo *> snpVec;

PatternSet patternUniqueTable;

bool verbose = false; // control phase reporting, statistics

// Progress messages.
void beginPhase()
{
  if (verbose)
  {
    //cout << "[" << msg << "... " << flush;
    time(&start_time);
  }
}

void beginPhase(const char *msg)
{
  if (verbose)
  {
    std::cout << "[" << msg << "... " << getpid() << std::endl;
    time(&start_time);
  }
}

void endPhase()
{
  if (verbose)
  {
    time_t end_time;
    time(&end_time);
    std::cout << " (" << difftime(end_time, start_time) << " s)]" << std::endl;
  }
}

void endPhase(const char *msg, std::string chr)
{
  if (verbose)
  {
    time_t end_time;
    time(&end_time);
    std::cout << "[" << msg << "... " << chr;
    std::cout << " (" << difftime(end_time, start_time) << " s)]" << std::endl;
  }
}

// show a pattern of any length.  Useful for the "merge" array used in combinePatterns
std::ostream &showPattern(std::ostream &os, const char *pattern, size_t len)
{
  for (int i = 0; i < (int)len; i++)
  {
    char c = pattern[i];
    if ('?' == c)
    {
      os << c;
    }
    else
    {
      os << (char)(c + '0');
    }
  }
  // for debugging
  os << std::flush;
  return os;
}

// show a pattern of length numStrains
std::ostream &showPattern(std::ostream &os, const char *pattern)
{
  return showPattern(os, pattern, numStrains);
}

// For debug print.
std::ostream &showPattern(const char *pattern)
{
  return showPattern(std::cout, pattern);
}

std::ostream &showPattern(const char *pattern, size_t len)
{
  return showPattern(std::cout, pattern, len);
}

// FIXME: This relevantStrains & strainAbbrevs should be "out" parameters, not globals.
// read in modified strain_index file with two columns: names (as before)
// and abbreviated names (for HTML).
// every strain should have an abbrev, whether we're working on it or not.
// But this was expedient.
void readStrains(char *fname)
{
  ColumnReader rdr(fname, (char *)"\t");

  int numtoks;
  while ((numtoks = rdr.getLine()) >= 0)
  {
    if (rdr.getCurrentLineNum() < 1) continue;
    // file has "Abbrev\tOfficalName\n"
    // if (numtoks != 2)
    // {
    //   std::cerr << "Warning: numtoks = " << numtoks << std::endl;
    // }

    // FIXME: some unnecessary string copies
    std::string name = rdr.getToken(1);
    std::string abbrev = rdr.getToken(0);

    relevantStrains.addElementIfNew(name);
    strainAbbrevs.addElementIfNew(abbrev);
  }
}


char *dupPattern(char *pattern)
{
    char *newpat = (char *)malloc(numStrains);
    std::memcpy(newpat, pattern, numStrains);
    return newpat;
}

// default constructor
SNPInfo::SNPInfo() : alleles(initAlleles()), frozen(false), used(false), qMarks(true) { pattern = initPattern(); };

// constructor for stuff in chr_info_perl file.
//  SNPInfo(string n, int ci, int pos) : name(n), chrIdx(ci), position(pos) { alleles = initAlleles(); };
SNPInfo::SNPInfo(std::string n, int ci, int pos) : name(n), chrIdx(ci), position(pos),
                                              alleles(initAlleles()),
                                              frozen(false), used(false), qMarks(true)
{
    pattern = initPattern();
};

// copy constructor
SNPInfo::SNPInfo(SNPInfo const &snpInfo) : name(snpInfo.name), chrIdx(snpInfo.chrIdx),
                                           position(snpInfo.position), alleles(strdup(snpInfo.alleles)),
                                           geneCodonMap(snpInfo.geneCodonMap),
                                           frozen(snpInfo.frozen), used(false), qMarks(true)
{

    std::cerr << "Illegal call to SNPInfo copy constructor" << std::endl;
    exit(1);
//     if (frozen) {
//       cerr << "Copy constructor applied to frozen SNPInfo:" << this << endl;
//     }
//     cout << "Pattern copy: old " << hex << (size_t) snpInfo.pattern << ", new "
// 	 << (size_t) pattern << endl;
//     pattern = dupPattern(snpInfo.pattern);
}

// destructor
SNPInfo::~SNPInfo()
{
    free(alleles);
    if (!frozen)
    {
        free(pattern);
    }
}

char * SNPInfo::initAlleles()
{
    if (numStrains < 0)
    {
        std::cerr << "FATAL: Cannot allocate allele strings until strains_index.txt has been read." << std::endl;
        exit(1);
    }
    char *tmpAlleles = (char *)malloc(numStrains + 1);
    std::memset(tmpAlleles, '?', numStrains);
    tmpAlleles[numStrains] = '\000';
    return tmpAlleles;
};

// has to be called after initAlleles.
char * SNPInfo::initPattern()
{
    char *tmpPattern = (char *)malloc(numStrains);
    std::memset(tmpPattern, '?', numStrains);
    return tmpPattern;
}

char SNPInfo::getAllele(int i) { return alleles[i]; }
void SNPInfo::setAllele(int i, char a) { alleles[i] = a; }



SNPInfo & SNPInfo::operator=(const SNPInfo &si)
{
    std::cerr << "Illegal call to SNPInfo assignment operator." << std::endl;
    exit(1);
    if (this != &si)
    {
        if (frozen)
        {
            std::cerr << "Assignment of frozen SNPInfo:" << this << std::endl;
        }
        name = si.name;
        chrIdx = si.chrIdx;
        position = si.position;
        qMarks = si.qMarks;
        std::strcpy(alleles, si.alleles);
        // memcpy(pattern, si.pattern, numStrains);
    }
    return *this;
}


void HaploBlock::init()
{
    pattern = (char *)malloc(numStrains);
    std::memset(pattern, '?', numStrains);
}

// default constructor
HaploBlock::HaploBlock() : chosen(false), score(0.0) { init(); }

// copy constructor
//
HaploBlock::HaploBlock(HaploBlock const &hb) : start(hb.start), size(hb.size),
                                   chosen(hb.chosen), score(hb.score)
{
    pattern = dupPattern(hb.pattern);
}

// Overload assignment to memcpy, so we don't have to worry about string copying.
HaploBlock & HaploBlock::operator=(const HaploBlock &hb)
{
    if (this != &hb)
    {
        start = hb.start;
        size = hb.size;
        chosen = hb.chosen;
        score = hb.score;
        std::memcpy(pattern, hb.pattern, numStrains);
    }
    return *this;
}

// destructor
HaploBlock::~HaploBlock() { free(pattern); }


// print a SNPInfo
std::ostream &operator<<(std::ostream &os, const SNPInfo &s)
{
    os << s.name << ", chromosome = " << chromosomes.eltOf(s.chrIdx)
       << ", location = " << s.position
       << ", alleles = " << s.alleles
       << ", pattern = ";
    showPattern(os, s.pattern);
    for (std::map<std::string, std::string>::const_iterator gmit = s.geneCodonMap.begin(); gmit != s.geneCodonMap.end(); gmit++)
    {
        os << ", " << (*gmit).first << " (" << (*gmit).second << ")";
    }
    os << std::endl;
    return os;
};

// print a HaploBlock
std::ostream &operator<<(std::ostream &os, const HaploBlock &hb)
{
    SNPInfo *pSNPInfo = snpVec[hb.start];
    os << "SNP " << hb.start << ": " << pSNPInfo->name
       << ", chromosome " << chromosomes.eltOf(pSNPInfo->chrIdx)
       << ", size " << hb.size
       << ", chosen  " << hb.chosen
       << ", pattern ";
    showPattern(hb.pattern);
    std::cout << ", score " << hb.score
         << std::endl;
    return os;
};


size_t PatternHash::operator()(const char *s) const
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

void showBlocks(std::vector<HaploBlock *> &hbv)
{
    for (unsigned i = 0; i < hbv.size(); i++)
    {
        std::cout << *(hbv[i]) << std::endl;
    }
}

// Test patterns for exact equality
bool eqPatterns(const char *pat1, const char *pat2)
{
    //  cout << "Pattern eq ";
    //  showPattern(pat1);
    //  cout << " = ";
    //  showPattern(pat2);
    bool result = std::memcmp(pat1, pat2, numStrains) == 0;
    //  cout << " result = " << result << endl;
    return result;
}