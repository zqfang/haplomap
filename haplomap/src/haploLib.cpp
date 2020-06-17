// -*-c++-*-

// Various useful functions from nhaploblocks

#include "haploLib.h"
#include <unistd.h>

static time_t start_time; // For time printing.  Should be in an object.

int numStrains = -1;
Dynum<string> relevantStrains;
Dynum<string> strainAbbrevs;

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
    cout << "[" << msg << "... " << getpid() << endl;
    time(&start_time);
  }
}

void endPhase()
{
  if (verbose)
  {
    time_t end_time;
    time(&end_time);
    cout << " (" << difftime(end_time, start_time) << " s)]" << endl;
  }
}

void endPhase(const char *msg, string chr)
{
  if (verbose)
  {
    time_t end_time;
    time(&end_time);
    cout << "[" << msg << "... " << chr;
    cout << " (" << difftime(end_time, start_time) << " s)]" << endl;
  }
}

// show a pattern of any length.  Useful for the "merge" array used in combinePatterns
ostream &showPattern(ostream &os, const char *pattern, size_t len)
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
  os << flush;
  return os;
}

// show a pattern of length numStrains
ostream &showPattern(ostream &os, const char *pattern)
{
  return showPattern(os, pattern, numStrains);
}

// For debug print.
ostream &showPattern(const char *pattern)
{
  return showPattern(cout, pattern);
}

ostream &showPattern(const char *pattern, size_t len)
{
  return showPattern(cout, pattern, len);
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
    // file has "Name\tAbbrev\n"
    if (numtoks != 2)
    {
      cerr << "Warning: numtoks = " << numtoks << endl;
    }

    // FIXME: some unnecessary string copies
    string name = rdr.getToken(1);
    string abbrev = rdr.getToken(0);

    relevantStrains.addElementIfNew(name);
    strainAbbrevs.addElementIfNew(abbrev);
  }
}


template <typename T>
decltype(std::bind(&T::value_type::second, std::placeholders::_1)) select2nd() {
    return std::bind(&T::value_type::second, std::placeholders::_1);
}


char *dupPattern(char *pattern)
{
    char *newpat = (char *)malloc(numStrains);
    memcpy(newpat, pattern, numStrains);
    return newpat;
}

// default constructor
SNPInfo::SNPInfo() : alleles(initAlleles()), frozen(false), used(false), qMarks(true) { pattern = initPattern(); };

// constructor for stuff in chr_info_perl file.
//  SNPInfo(string n, int ci, int pos) : name(n), chrIdx(ci), position(pos) { alleles = initAlleles(); };
SNPInfo::SNPInfo(string n, int ci, int pos) : name(n), chrIdx(ci), position(pos),
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

    cerr << "Illegal call to SNPInfo copy constructor" << endl;
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
        cout << "FATAL: Cannot allocate allele strings until strains_index.txt has been read." << endl;
        exit(1);
    }
    char *tmpAlleles = (char *)malloc(numStrains + 1);
    memset(tmpAlleles, '?', numStrains);
    tmpAlleles[numStrains] = '\000';
    return tmpAlleles;
};

// has to be called after initAlleles.
char * SNPInfo::initPattern()
{
    char *tmpPattern = (char *)malloc(numStrains);
    memset(tmpPattern, '?', numStrains);
    return tmpPattern;
}

SNPInfo & SNPInfo::operator=(const SNPInfo &si)
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
}


void HaploBlock::init()
{
    pattern = (char *)malloc(numStrains);
    memset(pattern, '?', numStrains);
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
        memcpy(pattern, hb.pattern, numStrains);
    }
    return *this;
}

// destructor
HaploBlock::~HaploBlock() { free(pattern); }


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

void showBlocks(vector<HaploBlock *> &hbv)
{
    for (unsigned i = 0; i < hbv.size(); i++)
    {
        cout << *(hbv[i]) << endl;
    }
}