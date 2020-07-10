//
// Created by Zhuoqing Fang on 7/8/20.
//

#ifndef HBCGM_EIGEN_H
#define HBCGM_EIGEN_H
#include <vector>
#include <string>
#include "dynum.h"

struct EigenMat {
    size_t size1;
    size_t size2;
    Dynum<std::string> rownames;
    gsl_matrix* data = nullptr; // original
    gsl_matrix* eigenvectors = nullptr;
    gsl_vector* eigenvalues = nullptr;
    gsl_vector* variances = nullptr;
    /// read correlation matrix(grm.rel) and rowname (grm.rel.id) files from PLink output,
    EigenMat(const char* MatrixFile, const char* RowNamesFile);
    /// init with STL vector
    explicit EigenMat(const std::vector<std::vector<double>> &Mat);
    /// init with GSL matrix
    explicit EigenMat(const gsl_matrix* Mat);
    EigenMat() = default;
    ~EigenMat();
    /// calculate eigenvectors and eigenvlaues
    void eigen();
    /// call eigen() first, and set negative value to 0
    void calcVariance();

    /// principal compoent analysis, return first L PCs
    gsl_matrix* pca(unsigned int L);
};

#endif //HBCGM_EIGEN_H