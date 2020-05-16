// -*-c++-*-

// (c) 2008, David L. Dill.  All rights reserved.

#ifndef UTIL_H
#define UTIL_H
#include <string>
#include <iostream>
// Random functions that should exist.

using namespace std;
inline const int toInt(const std::string& s)
{
  return std::atoi(s.c_str());
}

inline const float toFloat(const std::string& s)
{
  return std::atof(s.c_str());
}

// print vector  [1.0, 2.0, 3.0]

// Don't know how to get this to compile.
// print vector of strings ["foo", "bar", "bletch" ]
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << "[";
  typename vector<T>::const_iterator vend = v.end();
  for (typename vector<T>::const_iterator vit = v.begin(); vit < vend; vit++) {
    os << *vit;
    if (vit+1 < vend) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}

#endif