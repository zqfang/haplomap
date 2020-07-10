// -*-c++-*-

#ifndef DYNUM_H
#define DYNUM_H

// "Dynamic enum" datastructure.

// FIXME: add a checker, move some to non-inline, frozen bit
// FIXME: could make this more strongly typed by assigning
//   a class to the integers.  Could overload printing to print the original elt?
// but maybe that would be confusing.

#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>


template <typename EltType>
class Dynum
{
    typedef std::unordered_map<EltType, int> _E2IMap;
    _E2IMap _elt_to_idx;
    std::vector<EltType> _idx_to_elt;
    int _numIndices; // Number of indices allocated.  Also, the next index value.

public:
    Dynum() : _numIndices(0){}; // constructor
    ~Dynum(){}; // destructor

    // If elt is already in enum, return old idx.
    // Else, allocate a new idx and return it.
    int addElementIfNew(EltType elt)
    {
        typename _E2IMap::iterator it = _elt_to_idx.find(elt);
        if (it != _elt_to_idx.end())
        {
            // return *it.second;
            return it->second;
        }
        else
        {
            _idx_to_elt.push_back(elt);
            _elt_to_idx[elt] = _numIndices++;
            return _numIndices - 1;
        }
    }

    // Constructor to build a dynum out of an ordered vector of elements.
    // elements must be unique.
    Dynum(std::vector<EltType> &eltvec)
    {
        _idx_to_elt = eltvec; // FIXME: get rid of this copy.
        _numIndices = eltvec.size();
        for (int i = 0; i < eltvec.size(); i++)
        {
            // FIXME: check for already defined.
            _elt_to_idx[eltvec[i]] = i;
        }
    }

    // Map element -> index, but returns -1 if not indexed.
    int hasIndex(EltType elt)
    {
        typename _E2IMap::iterator it = _elt_to_idx.find(elt);
        if (it == _elt_to_idx.end())
        {
            return -1;
        }
        // return *it;
        return (*it).second;
    }

    // Map element -> index.  Fatal error if not there.
    // return -1 if not there
    int indexOf(EltType elt)
    {
        int idx = hasIndex(elt);
        if (idx >= 0)
        {
            return idx;
        }
        else
        {
            std::cerr << "Index of " << elt << " not found." << std::endl;
            return -1;
        }
    }

    // Map index -> element
    EltType eltOf(int idx)
    {
        return _idx_to_elt[idx];
    }

    // size
    int size()
    {
        return _numIndices;
    }

    // List mapping
    void dump()
    {
        for (int i = 0; i < size(); i++)
        {
            std::cout << i << "\t" << eltOf(i) << std::endl;
        }
    }
};

#endif
