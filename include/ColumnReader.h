// -*-c++-*-

// (c) 2008, David L. Dill.  All rights reserved.

// General-purpose function to read a file consisting of single-character-delimited columns.

#ifndef COLUMNREADER_H
#define COLUMNREADER_H


#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

using namespace std;

class ColumnReader {

  string _line;			// buffer for input line.
  vector<string> _line_vector; 	// vector of tokens for the current line.  Grows as necessary.
  char * _delimiters;			// column separator character.
  ifstream _in;			// input stream

  int _lineno;			// linenumber

public:

  ColumnReader(const char * fname, char *delimiters);		// constructor: open file.

  ~ColumnReader();		// destructor: close file.

  // reads line, fills _line_vector, returns number of tokens found.   
  // return -1 if no lines left in file.
  int getLine();		

  // Split a string at delimiter characters, fill token_vector with the substrings (clears token_vector).
  // General function -- perhaps it should go in another class.
  static int split(string s, char *delimiters, vector<string>& token_vector);

  string getToken(size_t i)
  {
    if (i < _line_vector.size()) {
      return _line_vector[i];
    }
    else {
      return string("");
    }
  }

  vector<string>::iterator begin() { return _line_vector.begin(); }; // return iterator to beginning of _line_vector

  vector<string>::iterator end() {return _line_vector.end(); } // return iterator to end of _line_vector

};

#endif

