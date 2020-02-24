// -*-c++-*-

// Various useful functions from nhaploblocks

#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>		// for numbers with commas in them.
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <iterator>
#include <algorithm>
#include <cmath>
//#include <ext/hash_set>
#include <ext/functional>
//#include <hash_map>

#include <time.h>

#include "util.h"
//#include "hash_string.h"
#include "dynum.h"
#include "indexcomparator.h"
#include "ColumnReader.h"

using namespace std;

// FIXME: tasteless global variables (maybe verbose is ok).

extern int numStrains;
extern bool verbose;		// in haploLib

extern Dynum<string> relevantStrains;
extern Dynum<string> strainAbbrevs;

ostream & showPattern(const char *pattern);
ostream & showPattern(const char *pattern, size_t len);
ostream & showPattern(ostream &os, const char *pattern);
ostream & showPattern(ostream &os, const char *pattern, size_t len);
void readStrains(char *fname);

// Progress messages.
void beginPhase();
void beginPhase(const char * msg);
void endPhase();
void endPhase(const char* msg, string chr);

