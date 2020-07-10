//
// Created by Zhuoqing Fang on 7/7/20.
//

#ifndef HBCGM_FDR_H
#define HBCGM_FDR_H

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <assert.h>
#include "ghmap.h"

/* Benjamini Hochberg procedure for controlling the FDR
 * sort pvalue in descending order, and return adjust pvalue in descending order.
 */
// see also: https://gist.github.com/wckdouglas/e3121c058c4fcf88cfd5
///Benjamini Hochberg procedure for controlling the FDR
///Takes a vector of p-values and alpha
///Return a vector of adjust-p-values
///Return a bool vector indicate p-values <= p are significant with FDR alpha.

template <typename _NumericType>
std::vector<bool> bh_fdr(std::vector<_NumericType>& pval, std::vector<_NumericType>& padj, float alpha)
{
    assert(pval.size() == padj.size());
    bool sorted = std::is_sorted(pval.begin(), pval.end(), std::greater<double>());
    /// decending order
    if (!sorted)
        std::sort(pval.begin(), pval.end(), std::greater<double>());

    // stored padj
    //std::vector<double> padj(pval.size(), 1.0);
    std::vector<bool> reject(pval.size(), false);
    _NumericType m = pval.size();
    uint32_t k = pval.size(); // This is the rank, doesn't need to be double.
    _NumericType factor;
    _NumericType p;
    for (int i=0; i < pval.size(); ++i) {
        factor = k/m;
        p = pval[i];
        if (p <= factor * alpha) {
            reject[i] = true;
        }
        p /= factor;
        padj[i] = p;
        k--; //Decrease rank
    }
    // accumulate minimun and make p = p < 1 ? p:1;
    std::accumulate(padj.begin(), padj.end(), 1.0, [](_NumericType x, _NumericType y){return std::min(x,y);});
    //std::for_each(padj.begin(),padj.end(),[](double &p){return p < 1 ? p: 1.0;});
    return reject;
}


/// flag: 1 sort pvalue, 0 sort mpvalue
std::vector<bool> bh_fdr(std::vector<BlockSummary *> & pval,
                         float alpha=0.05, bool flag = 1)
{
    std::vector<bool> reject(pval.size(), false);
    float m = pval.size();
    uint32_t k = pval.size(); // This is the rank, doesn't need to be double.
    float factor;
    float p;
    float previous_fdr;
    //BlockSummary* pBlock = new BlockSummary();
    // stored padj
    if (flag) {
        std::stable_sort(pval.begin(), pval.end(),
                  [](BlockSummary* x, BlockSummary* y) {return x->pvalue > y->pvalue;});
        previous_fdr =1.0;
        for (int i = 0; i < pval.size(); ++i) {
            factor = k / m;
            p = pval[i]->pvalue;
            if (p <= factor * alpha) {
                reject[i] = true;
            }
            p /= factor;
            pval[i]->FDR = std::min(p, previous_fdr); // accumulate minimum
            previous_fdr = p;
            k--; //Decrease rank
        }
//        std::accumulate(pval.begin(), pval.end(), pBlock,
//                        [](BlockSummary* x, BlockSummary* y)
//                        { return std::min(x->FDR, y->FDR); });
    } else {
        std::stable_sort(pval.begin(), pval.end(),
                  [](BlockSummary* x, BlockSummary* y) {return x->mPvalue > y->mPvalue;});
        previous_fdr = 1.0;
        for (int i = 0; i < pval.size(); ++i) {
            factor = k / m;
            p = pval[i]->mPvalue;
            if (p <= factor * alpha) {
                reject[i] = true;
            }
            p /= factor;
            pval[i]->mFDR = std::min(p, previous_fdr);
            previous_fdr = p;
            k--; //Decrease rank
        }
//        std::accumulate(pval.begin(), pval.end(), pBlock->mFDR,
//                        [](BlockSummary* x, BlockSummary* y)
//                                   { return std::min(x->mFDR, y->mFDR) ; });

    }
    //delete pBlock;
    return reject;
}








#endif //HBCGM_FDR_H
