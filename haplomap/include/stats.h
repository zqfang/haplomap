//
// Created by Zhuoqing Fang on 7/6/20.
//

#ifndef HBCGM_STATS_H
#define HBCGM_STATS_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cassert>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include "dynum.h"
#include "eigen.h"


/// ANONVA statistics
// This version takes vector of vectors in phenotypes, so that it can handle normal
// and categorical data in the same code.
// Handling of categorical data is as described to me by Ming.
// Means are means of vectors.
// SSW is sum of squares of Euclidean distances of phenotype vectors to the phenotype mean.
// SSB is sum of squares of Euclidean distances of haplotype averages to global mean.
// This does not actually compute the p-value.  It stops with the F statistic.
class ANOVA
{
private:
    std::vector<std::vector<float>> phenvec;
    unsigned int _numStrains;
    unsigned int _numHaplo;
    char * _pattern;
    float _FStat;
    float _pvalue; 
    float _effect;
    
    // sum vector, keep dimension
    std::vector<float> sumVector(std::vector<float> &vec);
    // subtract by a scalar, inplace
    void subVector(std::vector<float> &vec, float value);
    // element-wise addition inplace, destroys first argument (like +=)
    // v1 += v2
    void addVectors(std::vector<float> &v1, std::vector<float> &v2);
    // element-wise subtract inplace, v1 -= v2
    void subtractVectors(std::vector<float> &v1, std::vector<float> &v2);
    // element-wise multiply inplace, v1 *= v2
    float dotVectors(std::vector<float> &v1, std::vector<float> &v2);
    // multiply by scalar inplace.  Destroys first argument.
    void scaleVector(std::vector<float> &v1, float c);
    int numHaplotypes(char *pattern);
    int numDefinedStrains(char *pattern);

public:
    ANOVA(std::vector<std::vector<float>> &phenvec);
    void stat(char *pattern, float &FStat, float &pvalue, float &effect);
    ~ANOVA();
};


/// MANOVA analysis for population structure
class MANOVA
{
private:
    std::shared_ptr<EigenMat> _CorMat;
    std::vector<std::string> _MatRowNames;
    std::shared_ptr<Dynum<std::string>> _haploStrainsAbbrevs;
    gsl_matrix* _Mat; // matrix without qmark
    //gsl_matrix * _SubMat; // just submatrix view of _Mat
    bool _useEigen;
    unsigned int _numStrains;
    unsigned int _numDefined;
    unsigned int _numHaplo;
    unsigned int L;
    char* _pattern;

public:
    /// init with STL vector, set dimension reduction to L
    MANOVA(std::vector<std::vector<double>> &pc, unsigned int L=4);

    /// read correlation matrix(grm.rel) and rowname (grm.rel.id) files from PLink output,
    /// set dimension reduction to L
    MANOVA(const char* MatFile, char *delimiters, unsigned int L=4);
    /// init with GSL matrix, set dimension reduction to L
    MANOVA(gsl_matrix* M, unsigned int L=4);
    MANOVA() = default;
    ~MANOVA();

    /// return -1: pattern has too manny '?'
    /// return 0: success
    int setNonQMarkMat(char* pattern, Dynum<std::string>& haploStrainAbbr);

    /// calculate non negative eigenvalues
    void setEigen(Dynum<std::string>& haploStrainAbbr);

    /// get FSTAT and Pvalue using pilla's trace
    void pillaiTrace(float & FStat, float &PValue);

    /// read and read binary GSL matrix file
    void readMat(const char* filename, gsl_matrix * Mat);
    void writeMat(const char* filename, gsl_matrix * Mat);

    /// return new pattern without '?'. order remained.
    char * removeQMark(char *pattern);
    /// print haplotype pattern
    friend std::ostream& operator<<(std::ostream &os, const MANOVA& aov);



    /// returns maximum eq. class + 1.
    int numHaplotypes(char *pattern);

private:
    /// init _Mat, return new pattern without '?'
    void extractNonQMarkMat(gsl_matrix* M, char*pattern);

    /// drop Mat and pattern where strains marked with ?
    /// return new pattern without '?'
    void colMean(gsl_matrix * M, gsl_vector * mean);

    /// rows means and size group by pattern
    void groupMean(const gsl_matrix * M, char* pattern, gsl_matrix * _haploMean, gsl_matrix* _haploSize);

    /// total sum of squares
    void colTotalSumSquares(gsl_matrix * M, gsl_vector * tss);
    /// math trace
    double trace(const gsl_matrix * M);
    /// elemebt wise multiply, c*M
    void scale(gsl_matrix* M, double c);
    /// element wise sqrt
    void sqrt(gsl_matrix* M);
    /// dot product AB, save to C.
    void matmul(const gsl_matrix *A, const gsl_matrix *B, gsl_matrix* C);

};



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

/* Benjamini Hochberg procedure for controlling the FDR
 * sort pvalue in descending order, and return adjust pvalue in descending order.
 * flag: 1 sort pvalue, 0 sort mpvalue
 */
void bh_fdr(std::vector<BlockSummary *> & pval, float alpha=0.05, bool flag = 1);

#endif //HBCGM_STATS_H
