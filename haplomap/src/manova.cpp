//
// Created by Zhuoqing Fang on 7/6/20.
//
#include <stdio.h>
#include <string>
#include <cmath>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_linalg.h"
#include "ColumnReader.h"
#include "ghmap.h"
#include "manova.h"


MANOVA::MANOVA(const char* MatFile, const char* RowNameFile, unsigned int L):
_Mat(nullptr), _useEigen(false),_numStrains(0),_numDefined(0),_numHaplo(0),_pattern(nullptr)
{
    _CorMat = std::make_shared<EigenMat>(MatFile, RowNameFile);
    _numStrains = _CorMat->size1;
    assert(L > 0 && L < _CorMat->size2);
    this->L = L;
    // FIXME: select first L dim
    // use submatrix_view
//    gsl_matrix_view SubMat = gsl_matrix_submatrix(_CorMat->data, 0, 0, _numStrains, L);
//    _SubMat = &SubMat.matrix;
}

MANOVA::MANOVA(std::vector<std::vector<double>> & pc, unsigned int L):
_Mat(nullptr), _useEigen(false),_numStrains(0),_numDefined(0),_numHaplo(0),_pattern(nullptr)
{
    _CorMat = std::make_shared<EigenMat>(pc);
    assert(L > 0 && L < _CorMat->size2);
    this->L = L;
    // FIXME: select first L dim
    // use submatrix_view
//    gsl_matrix_view SubMat = gsl_matrix_submatrix(_CorMat->data, 0, 0, _numStrains, L);
//    _SubMat = &SubMat.matrix;
}

MANOVA::MANOVA(gsl_matrix* M, unsigned int L):
_Mat(nullptr), _useEigen(false),_numStrains(0),_numDefined(0),_numHaplo(0),_pattern(nullptr)
{
    _CorMat = std::make_shared<EigenMat>(M);
    assert(L > 0 && L < _CorMat->size2);
    this->L = L;
    // FIXME: select first L dim
    // use submatrix_view
//    gsl_matrix_view SubMat = gsl_matrix_submatrix(_CorMat->data, 0, 0, _numStrains, L);
//    _SubMat = &SubMat.matrix;
}

MANOVA::~MANOVA()
{
    if (_Mat != nullptr)
        gsl_matrix_free(_Mat);
}


void MANOVA::writeMat(const char* filename, gsl_matrix *Mat)
{
    FILE * f = fopen(filename, "wb");
    gsl_matrix_fwrite(f, Mat);
    fclose (f);
}

void MANOVA::readMat(const char* filename, gsl_matrix *Mat)
{
    FILE * f = fopen (filename, "rb");
    gsl_matrix_fread(f, Mat);
    fclose (f);
}


int MANOVA::numHaplotypes(char *pattern)
{
    int numHap = -1;
    for (int str1 = 0; str1 < _numStrains; str1++)
    {
        char hap = pattern[str1];
        if (pattern[str1] != '?' && numHap < hap)
        {
            numHap = hap ;
        }
    }
    return numHap + 1;
}

char* MANOVA::removeQMark(char *pattern)
{
    /// MARK:: now newpat could be returned
    char * newpat = (char *)malloc(_numStrains);
    std::memcpy(newpat, pattern, _numStrains);

    /// MARK: pattern is unprintable => strlen() = 0
    _numDefined = _numStrains;
    // remove all '?'
    int i = 0;
    while (i < _numDefined) {
        char hap = newpat[i];
        if (hap != '?') {
            // move left 1 step
            std::memmove(newpat+i, newpat+i+1, _numDefined - i);
            _numDefined --;
        } else {
            i++;
        }
    }
    //// debug reduced pattern
//    for (int i=0; i < _numDefined; ++i)
//        std::cout << (char)(newpatt[i]+'0'); // ASCII -> char
//    std::cout<<std::endl;
    //_pattern = _pat; // local memory
    return newpat;
}

void MANOVA::setEigen()
{
    assert(_CorMat != nullptr);
    _useEigen = true;
    _CorMat->eigen();
    _CorMat->calcVariance();

    /// explained variance
//    unsigned int var = 0;
//    for (int i = 0; i < _CorMat->size1; ++i)
//        if (5.0 < gsl_vector_get(_CorMat->variances, i))
//            var++;
//
//    this->L = std::min(var, this->L);

}

int MANOVA::setNonQMarkMat(char* pattern, Dynum<std::string>& haploStrainAbbr)
{
    //_numStrains = strlen(pattern);
    _numStrains = haploStrainAbbr.size();
    _haploStrainsAbbrevs = std::make_shared<Dynum<std::string>>(haploStrainAbbr);
    // called strdup to make makeUnprintable work
    _numHaplo = this->numHaplotypes(pattern);
    _pattern = pattern;

    // advoid memory leak
    if (_Mat != nullptr)
    {
        gsl_matrix_free(_Mat);
        _MatRowNames.clear();
        _Mat = nullptr;
        _numDefined = 0;
    }

    // define _numDefined, and get new pattern
    //this->removeQMark(pattern);
    for (int i = 0; i < _numStrains; i++)
    {
        if (pattern[i] != '?')
            _numDefined ++;
    }
    
    assert(_numDefined >= _numHaplo);
    // Stop run if
    if (_numDefined < _numStrains/2 || _numDefined < this->L)
    {
        std::cout<<"Skip too manny ? pattern: "<<this<<std::endl;
        return false;
    }

    //std::cout<<"residuals have less rank, MANOVA cannot be performed"<<std::endl;
    this->L = std::min(this->L, _numDefined - _numHaplo);

    gsl_matrix * pm;
    if (_useEigen)
    {
        assert(_CorMat->eigenvectors != nullptr);
        pm = _CorMat->eigenvectors;
    }
    else {
        pm = _CorMat->data;
    }
    this->extractNonQMarkMat(pm, pattern);

    return true;
}
void MANOVA::extractNonQMarkMat(gsl_matrix* M, char* pattern)
{
    // get new _Mat without '?'

    // re-assign
    //_Mat = gsl_matrix_alloc(_SubMat->size1, _SubMat->size2);
    _Mat = gsl_matrix_alloc(_numDefined, this->L);
    // skip '?' strains
    int rindex = 0, idx;
    std::string strain_abbr;
    for (int str1 = 0; str1 < _numStrains; str1++)
    {
        char hap = pattern[str1]; // 0,1,2,3,4, ?
        if ('?' != hap) {
            strain_abbr = _haploStrainsAbbrevs->eltOf(str1);
            //strains abbrevs exclude '?' ones
            _MatRowNames.push_back(strain_abbr);
            idx = _CorMat->rownames.indexOf(strain_abbr);
            for (int cindex=0; cindex < this->L; ++cindex)
            {
                double value = gsl_matrix_get(M, idx, cindex);
                gsl_matrix_set(_Mat, rindex, cindex, value);
            }
            rindex ++;
        }
    }
    assert(_numDefined == rindex);

}

void MANOVA::pillaiTrace(float & FStat, float &PValue )
{
    assert(_Mat != nullptr);
    int cols = _Mat->size2;
    int rows = _Mat->size1;

    // Compute means groupby haplotype
    gsl_matrix* haploMean = gsl_matrix_alloc(_numHaplo, cols);
    gsl_matrix_set_zero(haploMean);
    gsl_matrix* haploSize = gsl_matrix_alloc(_numHaplo, cols);
    gsl_matrix_set_zero(haploSize);
    groupMean(_Mat, _pattern, haploMean, haploSize);

    // make a copy
    // gsl_matrix* _data = gsl_matrix_alloc(rows, cols);
    // gsl_matrix_memcpy(_data, _Mat);
    // gsl_matrix_free(_data);

    // dimension not match!
    // gsl_vector* tss = gsl_vector_alloc(cols);
    // colTotalSumSquares(_data, tss);

    // grand mean
    gsl_vector* grandMean = gsl_vector_alloc(cols);
    colMean(_Mat, grandMean);

    gsl_matrix* _grandMean0 = gsl_matrix_alloc(rows, cols);
    for (int i = 0; i< rows; i++)
        gsl_matrix_set_row(_grandMean0, i, grandMean);
    // Y_{ij} - Y_bar
    gsl_matrix_sub(_Mat, _grandMean0);
    gsl_matrix_free(_grandMean0);
    // total sum of squares
    gsl_matrix* tss = gsl_matrix_alloc(_Mat->size2, _Mat->size2);
    // matmul
    gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                   1.0, _Mat, _Mat, 0.0, tss);

    gsl_matrix* _grandMean = gsl_matrix_alloc(haploMean->size1, haploMean->size2);
    for (int i = 0; i< _numHaplo; i++)
        gsl_matrix_set_row(_grandMean, i, grandMean);
    // Y_{j} - Y_bar
    gsl_matrix_sub(haploMean, _grandMean);
    gsl_matrix_free(_grandMean);
    // sqrt(haploSize)
    this->sqrt(haploSize);
    gsl_matrix_mul_elements(haploMean, haploSize);

    // hypothesis sum of squares matrix
    gsl_matrix* hss = gsl_matrix_alloc(haploMean->size2, haploMean->size2);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                   1.0, haploMean, haploMean, 0.0, hss);

    // ess = tss - hss
    // Pillai’s Trace
    // trace tr(H / (E+H)) =  sum(diag(hss/tss))
    gsl_matrix_div_elements(hss, tss);
    double vv = this->trace(hss);

    //int df_t = rows -1;
    int df_h = _numHaplo -1;
    int e = rows - _numHaplo;
    int p = cols;
    int s = std::min(p, df_h);
    double m = (std::abs(p-df_h) - 1 ) /2.0;
    // FIXME: e > p+1?
    //assert(e >= p);
    if (e < p) {
        FStat = INFINITY;
        PValue = 1.0;

        // debugging
        // could not process
        std::cerr<<"numHaplo: "<<_numHaplo
                 <<" numDefine: "<<_numDefined
                 <<" numStrains: "<<_numStrains
                 <<" pattern: ";

        for (int i=0; i < _numStrains; ++i) {
            char hap = _pattern[i];
            if (hap != '?')
                std::cerr << (char) (hap + '0'); // ASCII -> char
            else
                std::cerr << hap;
        }
        std::cerr<<std::endl;
        return;
    }
    double n = (e - p - 1) /2.0;
    // FIXME:: penalize haplotypes with ?
    //1.1*(_numStrains - _numDefined);
    double df_num = s*(2*n +s +1);
    double df_den = s*(2*m+s+1);
    double f_num = (2*n+s+1)*vv;
    double f_den = (2*m+s+1)*(s-vv);
    FStat = (float) ( f_num / f_den);
    PValue = (float) gsl_cdf_fdist_Q(FStat, df_den, df_num);
    if (gsl_isnan(PValue))
        PValue =1.0;

    gsl_vector_free(grandMean);
    gsl_matrix_free(haploMean);
    gsl_matrix_free(haploSize);
    gsl_matrix_free(hss);
    gsl_matrix_free(tss);
}

void MANOVA::sqrt(gsl_matrix * M) {
    for(int i = 0; i < M->size1; ++i)
    {
        for (int j = 0; j< M->size2; ++j)
        {
            double x = gsl_matrix_get(M, i, j);
            gsl_matrix_set(M, i,j,std::sqrt(x));
        }
    }
}

double MANOVA::trace(const gsl_matrix *M) {
    double vv = 0;
    for (int i = 0; i< M->size1;++i)
        vv += gsl_matrix_get(M,i,i);

    return vv;
}

void MANOVA::groupMean(const gsl_matrix *M, char *pattern, gsl_matrix* _haploMean, gsl_matrix* _haploSize) {

    gsl_matrix_set_zero(_haploMean);
    gsl_matrix_set_zero(_haploSize);

    // FIXME: assert(M->size1 == strlen(pattern));
    for (int j = 0; j < M->size2; j++)
    {
        for (int i = 0; i < M->size1; i++)
        {
            char hap = pattern[i]; // 0,1,2,3,4, ?
            if ('?' != hap)
            {
                // FIXME: char to int ?
                double m = gsl_matrix_get(_haploMean, hap, j);
                m += gsl_matrix_get(M, i, j);
                gsl_matrix_set(_haploMean, hap, j, m);
                //
                double n = gsl_matrix_get(_haploSize, hap, j);
                n ++;
                gsl_matrix_set(_haploSize, hap, j, n);
            }
        }
    }
    // get real mean
    gsl_matrix_div_elements(_haploMean, _haploSize);
}

void MANOVA::colMean(gsl_matrix * M, gsl_vector * mean)
{
    gsl_vector_set_zero(mean);
    // calculate the mean of each column
    for (int j=0; j < M->size2; j++)
    {
        //get current column as vector
        gsl_vector_view myColumn = gsl_matrix_column(M, j);
        // calc mean
        double avg = gsl_stats_mean(myColumn.vector.data,myColumn.vector.stride, myColumn.vector.size);
        // assign
        gsl_vector_set(mean, j, avg);
    }
}

void MANOVA::colTotalSumSquares(gsl_matrix *M, gsl_vector *tss) {
    gsl_vector_set_zero(tss);
    // calculate the mean of each column
    for (int j=0; j < M->size2; j++) {
        //get current column as vector
        gsl_vector_view myColumn = gsl_matrix_column(M, j);
        double t = gsl_stats_tss(myColumn.vector.data,myColumn.vector.stride, myColumn.vector.size);
        // assign
        gsl_vector_set(tss, j, t);
    }
}

void MANOVA::scale(gsl_matrix *M, double c) {
    assert(M != NULL);
    unsigned int rows = M->size1;
    unsigned int cols = M->size2;

    for (int i=0; i < rows; i++){
      for(int j=0; j < cols; j++){
          double myDobule = gsl_matrix_get(M, i, j);
          myDobule *= c;
          // update matrix
          gsl_matrix_set(M, i, j, myDobule);
      }
    }
}

// matrix multiplication GSL_BLAS
// https://www.gnu.org/software/gsl/doc/html/blas.html
// s,d,c, z = double precision ...
// ge = general matrices
// mm = matrix-matrix multiplication
void MANOVA::matmul(const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *C)
{
    // OUT = alpha* op(A)op(B) + beta*C
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, A, B,
                   0.0, C);
}

std::ostream& operator<<(std::ostream &os, const MANOVA& aov)
{
    // for debugging
    for (int i=0; i < aov._numStrains; ++i)
    {
        char hap = aov._pattern[i];
        if (hap != '?')
            os << (char)(hap + '0'); // ASCII -> char
        else
            os << hap;
    }
    return os;
}


namespace HBCGM {
    // debugging
    void print_gvec(gsl_vector* M)
    {
        std::cout<<"trace vector: "<<std::endl;
        for (int i = 0; i < M->size; ++i)
        {
            std::cout<<gsl_vector_get(M, i)<<" ";
        }
        std::cout<<std::endl;
    }

    void print_gmat(gsl_matrix* M){
        std::cout << "trace matrix: "<<std::endl;
        for (int i = 0; i < M->size1; ++i){
            for (int j=0; j < M->size2; ++j) {
                std::cout<< gsl_matrix_get(M, i, j) << " ";
            }
            std::cout<<std::endl;
        }
    }
}
