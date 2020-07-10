// -*-c++-*-

#ifndef INDEXCOMPARATOR_H
#define INDEXCOMPARATOR_H

#include <fstream>
#include <ostream>
#include <vector>
#include <functional>

//using namespace std;

// Function object for sorting an array of indices based on another array of values.

// *** Need to make this into a template

template <typename sortType, typename lss>
struct IndexComparator : public std::binary_function<sortType, sortType, bool>
{
  std::vector<sortType> *_pCmpVec;

  IndexComparator(std::vector<sortType> *pcmpVec) : _pCmpVec(pcmpVec){};

  bool operator()(const int &a, const int &b)
  {
    //    cout << "Compare a = " << a << "(" << (*_pCmpVec)[a] << ") vs. ";
    //    cout << "b = " << b << "(" << (*_pCmpVec)[b] << ")" << " -> ";
    // bool result = ((*_pCmpVec)[a] > (*_pCmpVec)[b]);
    bool result = lss()((*_pCmpVec)[a], (*_pCmpVec)[b]);
    //    cout << result << endl;
    return result;
  };
};

// Usage:

//   // Sort a strOrderVec of indices 0..numStrains-1 so that they are in descending order
//   // of numDefVec.  Color the leftmost strains first.

//   vector<int>::iterator stoEnd = strOrderVec.end();
//   int i = 0;
//   for (vector<int>::iterator stoIt = strOrderVec.begin(); stoIt != stoEnd; stoIt++) {
//     *stoIt = i++;
//   }

//   IndexComparator idxCompare(&numDefVec);

//   stable_sort(strOrderVec.begin(), strOrderVec.end(), idxCompare);

#else
#endif
