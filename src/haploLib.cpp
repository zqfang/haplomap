// -*-c++-*-

// Various useful functions from nhaploblocks

#include "haploLib.h"
#include <sys/types.h>
#include <unistd.h>

static time_t start_time; // For time printing.  Should be in an object.

int numStrains = -1;
Dynum<string> relevantStrains;
Dynum<string> strainAbbrevs;

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
    string name = rdr.getToken(0);
    string abbrev = rdr.getToken(1);

    relevantStrains.addElementIfNew(name);
    strainAbbrevs.addElementIfNew(abbrev);
  }
}
