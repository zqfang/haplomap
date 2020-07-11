//
// Created by Zhuoqing Fang on 7/6/20.
//

#ifndef HBCGM_MANOVA_H
#define HBCGM_MANOVA_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include "dynum.h"
#include "eigen.h"


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
    MANOVA(const char* MatFile, const char* RowNameFile, unsigned int L=4);

    /// init with GSL matrix, set dimension reduction to L
    MANOVA(gsl_matrix* M, unsigned int L=4);
    MANOVA() = default;
    ~MANOVA();

    /// return -1: pattern has too manny '?'
    /// return 0: success
    int setNonQMarkMat(char* pattern, Dynum<std::string>& haploStrainAbbr);

    /// calculate non negative eigenvalues
    void setEigen();

    /// get FSTAT and Pvalue using pilla's trace
    void pillaiTrace(float & FStat, float &PValue);

    /// read and read binary GSL matrix file
    void readMat(const char* filename, gsl_matrix * Mat);
    void writeMat(const char* filename, gsl_matrix * Mat);

    /// print haplotype pattern
    friend std::ostream& operator<<(std::ostream &os, const MANOVA& aov);


private:
    /// returns maximum eq. class + 1.
    int numHaplotypes(char *pattern);
    /// strip all '?'
    void removeQMark(char *pattern);

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
#endif //HBCGM_MANOVA_H
