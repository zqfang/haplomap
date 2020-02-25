// -*-c++-*-

// (c) 2008, David L. Dill.  All rights reserved.

// General-purpose function to read a file consisting of single-character-delimited columns.

#include "../include/ColumnReader.h"

ColumnReader::ColumnReader(const char * fname, char *delimiters) {
  _in.open(fname);
  _delimiters = delimiters;
  if (!_in.is_open()) {
    cerr << "Open of file \"" << fname << "\" failed: ";
    perror("");
    exit(1);
  }
  _lineno = 0;
};

ColumnReader::~ColumnReader() {
  _in.close();
}



// General function -- perhaps it should go in another class.
// Similar to perl split
// int ColumnReader::split(string& s, char *delimiters, vector<string>& token_vector)
int ColumnReader::split(string s, char *delimiters, vector<string>& token_vector)
{
  token_vector.clear();

  size_t tokstart = 0;
  size_t dpos = 0;

  // boundary cases:
  // if s is empty, loop does nothing, vector will be empty.
  // if s has one thing, ffo returns npos, first thing gets pushed,
  // then it returns.
  // This was segfaulting; I'm wondering if that's a bug in substr.
  // In 64-bit mode, ffo seems to return a very large number.  But npos
  // appears to be a 32-bit quantity.
  while (dpos < s.size()) {
    dpos = s.find_first_of(delimiters, tokstart);
    // this crock necessitated by weird segfault.
    size_t sz = s.size();
    if (dpos > sz) {
      dpos = sz;
    }
    string tok = s.substr(tokstart, dpos-tokstart);
    token_vector.push_back(tok);
    tokstart = dpos+1;		// don't care if it overflows
  }
  
  return token_vector.size();
}

int ColumnReader::getLine() 
{
  _lineno++;
  if (getline(_in, _line).eof()) {
    // no lines remain.
    return -1;
  }

  // tokenize it
  return split(_line, _delimiters, _line_vector);
}

