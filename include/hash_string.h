// -*-c++-*-

#ifndef HASH_STRING_H
#define HASH_STRING_H

// (c) 2008, David L. Dill.  All rights reserved.

// Annoyingly, the gnu C++ library does not supply a standard hash function for strings.
// So this is it.

//#include <ext/hash_map>

namespace __gnu_cxx {
template<> struct hash< std::string >                                                       
{                                                                                           
  size_t operator()( const std::string& x ) const                                           
  {                                                                                         
    return hash< const char* >()( x.c_str() );                                              
  }                                                                                         
};                                                                                          
}

#else
#endif

