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



class ColumnReader
{
    std::string _line;                // buffer for input line.
    std::vector<std::vector<std::string>> _header_vectors;  // header lines, store multi-header-lines which start with "#"
    std::vector<std::string> _line_vector; // vector of tokens for the current line.  Grows as necessary.
    char *_delimiters;           // column separator character.
    std::ifstream _in;                // input stream
    int _lineno; // linenumber indicator. Note:: header lines will maintain _lineno = 0

public:
    ColumnReader(const char *fname, char *delimiters); // constructor: open file.
    ~ColumnReader(); // destructor: close file.
    // reads line, fills _line_vector, returns number of tokens found.
    // return -1 if no lines left in file.
    int getLine();
    int getCurrentLineNum();
    std::string getCurrentLineString();
    std::vector<std::vector<std::string>> getHeaderLines(void);
    // Split a string at delimiter characters, fill token_vector with the substrings (clears token_vector).
    // General function -- perhaps it should go in another class.
    static int split(std::string s, char *delimiters, std::vector<std::string> &token_vector);
    std::string strip(const std::string& str, const std::string delimiter = " \n\r\t");
    std::string join(std::vector<std::string> &, char * delimiter);
    std::string getToken(size_t i);
    std::vector<std::string>::iterator begin() { return _line_vector.begin(); }; // return iterator to beginning of _line_vector
    std::vector<std::string>::iterator end() { return _line_vector.end(); } // return iterator to end of _line_vector
};

#endif
