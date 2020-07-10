//
// Created by Zhuoqing Fang on 7/8/20.
//

#ifndef HBCGM_ARGSORT_H
#define HBCGM_ARGSORT_H
#include <vector>


/// argsort usage:
///for (auto i: argsort(v))
///{
///   cout << v[i] << endl;
///}
template <typename T>
std::vector<size_t> argsort(const std::vector<T> &v)
{

// initialize original index locations
    std::vector<size_t> idx(v.size());
// Fills the range [first, last) with sequentially increasing values,
// starting with value and repetitively evaluating ++value
    std::iota(idx.begin(), idx.end(), 0);

// sort indexes based on comparing values in v
// using std::stable_sort instead of std::sort
// to avoid unnecessary index re-orderings
// when v contains elements of equal values
    std::stable_sort(idx.begin(), idx.end(),
                     [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}


template <typename T, typename _Compare>
std::vector<size_t> argsort(const std::vector<T> &v, _Compare comp)
{

// initialize original index locations
    std::vector<size_t> idx(v.size());
// Fills the range [first, last) with sequentially increasing values,
// starting with value and repetitively evaluating ++value
    std::iota(idx.begin(), idx.end(), 0);

// sort indexes based on comparing values in v
// using std::stable_sort instead of std::sort
// to avoid unnecessary index re-orderings
// when v contains elements of equal values
    std::stable_sort(idx.begin(), idx.end(),
                     [&v,&comp](size_t i1, size_t i2) -> bool {return comp(v[i1],v[i2]);});
    return idx;
}


template <class _RAIter, class _Compare>
std::vector<size_t> argsort(_RAIter first, _RAIter last, _Compare comp) {

    std::vector<size_t> idx(last-first);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
                     [&first,comp](size_t i1, size_t i2) {return comp(first[i1], first[i2]);});

    return idx;
}



#endif //HBCGM_ARGSORT_H
