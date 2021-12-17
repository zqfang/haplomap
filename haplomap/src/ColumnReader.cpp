// -*-c++-*-

// (c) 2008, David L. Dill.  All rights reserved.

// General-purpose function to read a file consisting of single-character-delimited columns.

#include "ColumnReader.h"


ColumnReader::ColumnReader(const char *fname, char *delimiters):
_delimiters(delimiters)
{
    //_in.imbue(std::locale());
    _in.open(fname, std::ios::in);
    _delimiters = delimiters;
    if (!_in.is_open())
    {
      std::cerr << "Open of file \"" << fname << "\" failed: " << std::endl;
      std::exit(1);
    }
    _lineno = 0;

    // if (std::getline(_in, _line).eof())
    // {
    //   // no lines remain.
    //   std::cerr << "Empty Input file: "<<fname<<std::endl;
    //   std::exit(1);
    // }

}

ColumnReader::~ColumnReader()
{
    _in.close();
}

// General function -- perhaps it should go in another class.
// Similar to perl split
int ColumnReader::split(std::string s, char *delimiters, std::vector<std::string> &token_vector)
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
  while (dpos < s.size())
  {
    dpos = s.find_first_of(delimiters, tokstart);
    // this crock necessitated by weird segfault.
    size_t sz = s.size();
    if (dpos > sz)
    {
      dpos = sz;
    }
    std::string tok = s.substr(tokstart, dpos - tokstart);
    token_vector.push_back(tok);
    tokstart = dpos + 1; // don't care if it overflows
  }

  return token_vector.size();
}


/// trim lines
std::string ColumnReader::strip(const std::string& str, const std::string delimiter)
{
//    std::string s;
//    s.erase(s.find_last_not_of(" \n\r\t")+1);
   size_t first = str.find_first_not_of(delimiter);
   if (std::string::npos == first)
   {
       return str;
   }
   size_t last = str.find_last_not_of(delimiter);
   return str.substr(first, (last - first + 1));
}


std::vector<std::vector<std::string>> ColumnReader::getHeaderLines(void)
{
  return _header_vectors;
}

int ColumnReader::getCurrentLineNum()
{
  return _lineno;
}
std::string ColumnReader::getCurrentLineString()
{
  return _line;
}
int ColumnReader::getLine()
{
  if (getline(_in, _line).eof())
  {
    // no lines remain.
    return -1;
  }

  if (_line.find("#") == 0)
  {
      // comment/header line stats with #
      size_t start = _line.find_first_not_of('#');
      std::vector<std::string> __head;
      int num = this->split(_line.substr(start), _delimiters, __head);
      this->_header_vectors.push_back(__head);
      return num;
    }
  _lineno++;
  // tokenize it
  return this->split(_line, _delimiters, _line_vector);
}


std::string ColumnReader::getToken(size_t i) 
{
  if (i < _line_vector.size())
  {
  return _line_vector[i];
  }
  else
  {
  return std::string("");
}
}
