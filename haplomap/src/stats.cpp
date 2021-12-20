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
#include "gsl/gsl_math.h"
#include "ColumnReader.h"
#include "ghmap.h"
#include "stats.h"

void bh_fdr(std::vector<BlockSummary *> & pval, float alpha, bool flag)
{
    //std::vector<bool> reject(pval.size(), false);
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
        for (unsigned i = 0; i < pval.size(); ++i) {
            factor = k / m;
            p = pval[i]->pvalue;
            //if (p <= factor * alpha) {
            //    pval[i]->relReject = true;
            //}
            p /= factor;
            pval[i]->FDR = std::min(p, previous_fdr); // accumulate minimum
            previous_fdr = pval[i]->FDR;
            k--; //Decrease rank
        }
//        std::accumulate(pval.begin(), pval.end(), pBlock,
//                        [](BlockSummary* x, BlockSummary* y)
//                        { return std::min(x->FDR, y->FDR); });
    } else {
        std::stable_sort(pval.begin(), pval.end(),
                         [](BlockSummary* x, BlockSummary* y) {return x->relPvalue > y->relPvalue;});
        previous_fdr = 1.0;
        for (unsigned i = 0; i < pval.size(); ++i) {
            factor = k / m;
            p = pval[i]->relPvalue;
            if (p <= factor * alpha) {
                pval[i]->relReject = true;
            }
            p /= factor;
            pval[i]->relFDR = std::min(p, previous_fdr);
            previous_fdr = pval[i]->relFDR;
            k--; //Decrease rank
        }
//        std::accumulate(pval.begin(), pval.end(), pBlock->mFDR,
//                        [](BlockSummary* x, BlockSummary* y)
//                                   { return std::min(x->mFDR, y->mFDR) ; });

    }
    //delete pBlock;
}




ANOVA::ANOVA(std::vector<std::vector<float>> &phenvec):
_numHaplo(0),_pattern(nullptr),_FStat(INFINITY),_pvalue(1.0),_effect(0)
{
  this->phenvec = phenvec;
  this->_numStrains = phenvec.size();
}
ANOVA::~ANOVA(){}

std::vector<float> ANOVA::sumVector(std::vector<float> &vec) {
    std::vector<float> res(1,0.0F);
    //res[0] = std::accumulate(vec.begin(), vec.end(), 0.0F);
    for (float & v: vec )
        res[0] += v;
    return res;
}

void ANOVA::subVector(std::vector<float> &vec, float value) {
    for (float & v: vec ) 
        v -= value;
}

// Some vector arithmetic.
// destroys first argument (like +=)
void ANOVA::addVectors(std::vector<float> &v1, std::vector<float> &v2)
{
    if (v1.size() != v2.size())
    {
        std::cout << "addVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << std::endl;
        exit(1);
    }
    std::vector<float>::iterator vend = v1.end();
    std::vector<float>::iterator vit2 = v2.begin();
    for (std::vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
    {
        *vit1 += *vit2;
        vit2++;
    }
}

void ANOVA::subtractVectors(std::vector<float> &v1, std::vector<float> &v2)
{
    if (v1.size() != v2.size())
    {
        std::cout << "subtractVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() <<std::endl;
        exit(1);
    }
    std::vector<float>::iterator vend = v1.end();
    std::vector<float>::iterator vit2 = v2.begin();
    for (std::vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
    {
        *vit1 -= *vit2;
        vit2++;
    }
}

//
float ANOVA::dotVectors(std::vector<float> &v1, std::vector<float> &v2)
{
    float result = 0.0;
    if (v1.size() != v2.size())
    {
        std::cout << "dotVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << std::endl;
        exit(1);
    }
    std::vector<float>::iterator vend = v1.end();
    std::vector<float>::iterator vit2 = v2.begin();
    for (std::vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
    {
        result += (*vit1) * (*vit2);
        vit2++;
    }
    return result;
}

// multiply by scalar.  Destroys first argument.
void ANOVA::scaleVector(std::vector<float> &v1, float c)
{
    std::vector<float>::iterator vend = v1.end();
    for (std::vector<float>::iterator vit = v1.begin(); vit < vend; vit++)
    {
        *vit *= c;
    }
}
int ANOVA::numHaplotypes(char *pattern)
{
    int numHap = -1;
    // int _numStrains = strlen(pattern);
    for (unsigned str1 = 0; str1 < this->_numStrains; str1++)
    {
        int hap = pattern[str1];
        if (pattern[str1] != '?' && numHap < hap)
        {
            numHap = hap ;
        }
    }
    return numHap + 1;
}

int ANOVA::numDefinedStrains(char *pattern)    
{
    int count = 0;
    for (unsigned str1 = 0; str1 < this->_numStrains; str1++)
    {
        if (pattern[str1] != '?')
        {
            count++;
        }
    }
    return count;
}

void ANOVA::stat(char *pattern, float &FStat, float &pvalue, float &effect)
{
    int _numHaplo = this->numHaplotypes(pattern);
    // int _numStrains = strlen(pattern);
    // array haplotype -> num strains in haplotype.
    std::vector<int> haploNum(_numHaplo, 0);
    // array haplotype -> mean (std::vector<float>) for each haplotype. Note: numCategories is a global variable
    std::vector<std::vector<float>> haploMean(_numHaplo, std::vector<float>(numCategories, 0.0F)); // size 1

    float numDefined = 0.0F; // numbers of strains without ?
    float numDefinedDataPoints = 0.0F; // number of all individual data point without '?' strain
    std::vector<float> sumDefined(numCategories, 0.0F); //size 1

    std::vector<std::vector<float>> sumStrains;

    // Compute haplotype means
    for (unsigned str1 = 0; str1 < _numStrains; str1++)
    {
        int hap = pattern[str1]; // 0,1,2,3,4
        std::vector<float> &phen = this->phenvec[str1]; // multivalue?
        if ('?' != hap)
        {
            numDefined ++;
            numDefinedDataPoints += phen.size();
            sumStrains.push_back(sumVector(phen));
            //addVectors(sumDefined, sumStrains.back());
            //haploNum[hap]++;

            //addVectors(haploMean[hap], phen); // temporarily, the total, not mean.
            addVectors(haploMean[hap], sumStrains.back());
            haploNum[hap] += phen.size();
        }
    }

    for (int hap = 0; hap < _numHaplo; hap++)
    {
        scaleVector(haploMean[hap], 1.0 / haploNum[hap]); // get the mean
    }

    if (traceFStat)
    {
        std::cout << "numDefined = " << numDefined << ", sumDefined = " << sumDefined.back() << std::endl;
        std::cout << "haploNum[] = [";
        for (int hap = 0; hap < _numHaplo; hap++)
        {
            std::cout << haploNum[hap] << " ";
        }
        std::cout << "]" << std::endl;

        std::cout << "haploMean[] = [";
        for (int hap = 0; hap < _numHaplo; hap++)
        {
            std::cout << haploMean[hap].back() << " ";
        }
        std::cout << "]" << std::endl;
    }

    std::vector<float> mean(1, 0.0F); // mean of all data
    for (auto & s: sumStrains)
        mean[0] += s.back();
    // scaleVector(mean, 1.0 / numDefined);
    scaleVector(mean, 1.0 / numDefinedDataPoints);

    //float SST = 0.0F;
    float SSW = 0.0F;
    for (unsigned str1 = 0; str1 < _numStrains; str1++)
    {
        int hap = pattern[str1];
        if ('?' != hap)
        {
            std::vector<float> resid = phenvec[str1];
            subVector(resid, haploMean[hap].back());
            SSW += dotVectors(resid, resid);
            // vector<float> resid2 = phenvec[str1];
            // subVector(resid2, mean.back());
            // SST += dotVectors(resid2, resid2);
        }
    }

    // SSB -- between sum of squares (sum over haplotypes hapsize*(hapmean-mean)^2
    float SSB = 0.0;
    for (int hap = 0; hap < _numHaplo; hap++)
    {
        std::vector<float> diff = haploMean[hap]; // copy so we don't destroy haploMeans
        subtractVectors(diff, mean);         // (haplotype mean) - mean
        float sq = haploNum[hap] * dotVectors(diff, diff);
        SSB += sq;
    }
    // FIXME:
    // my quick and dirty hack to penalize missing alleles.
    // simulates additional error for each missing value (but ignores
    // degrees of freedom).
    SSW += 1.1 * (_numStrains - numDefined);

    // mean square within
    // df within is (numDefined-numHaplo)
    float dfW = numDefinedDataPoints - _numHaplo; // degrees of freedom within
    float dfB = _numHaplo - 1;          // degrees of freedom between.
    if (dfW == 0.0)
    {
        // This happens when numDefined = numHaplo, which occurs rarely when there are
        // lots of undefined strains and lots of haplotypes.
        // This causes errors in gsl, and I don't know what the right thing to do is,
        // so just punt.  We won't get a match for this.
        pvalue = 1.0;
        effect = 0.0;
        return;
    }
    float MSW = SSW / dfW;
    float MSB = SSB / dfB;

    // This formula is the same as Peltz.  So, I think I've seen two totally
    // different formulas for omega^2
    // WARNING: This divides by 0 if SSW is 0.
    // Which seems to work ok (F <- "inf").
    FStat = MSB / MSW;

    // out parameter for pvalue
    if (gsl_isnan(FStat)) {
        std::cout<<"FStat is NaN for entry: "<<" numHaplo: "
        << _numHaplo <<", numStrains: "<< _numStrains << ", numDefinedStrains: " << numDefined <<", Haplotype: ";
        for (unsigned i=0; i < _numStrains; ++i)
            std::cout << (char)(pattern[i]+'0'); // ASCII -> char
        std::cout<<std::endl;

    }
    pvalue = (float)gsl_cdf_fdist_Q((double)FStat,(double)dfB,(double)dfW);
    if (gsl_isnan(pvalue))
        pvalue = 1.0;

    // Genetic effect
    // Oh wow!  omega^2 is the "coefficient of determination"!
    // http://faculty.chass.ncsu.edu/garson/PA765/anova.htm#anova2
    effect = (float)((SSB - ( _numHaplo - 1) * MSW) / (SSW + SSB + MSW));
}



MANOVA::MANOVA(const char* MatFile, char* delimiter, unsigned int L):
_Mat(nullptr), _useEigen(false),_numStrains(0),_numDefined(0),_numHaplo(0),_pattern(nullptr)
{
    _CorMat = std::make_shared<EigenMat>(MatFile, delimiter);
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

    if (_pattern != nullptr)
        free(_pattern);
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
    for (unsigned str1 = 0; str1 < _numStrains; str1++)
    {
        int hap = pattern[str1];
        if (pattern[str1] != '?' && numHap < hap)
        {
            numHap = hap ;
        }
    }
    return numHap + 1;
}

char* MANOVA::removeQMark(char *pattern)
{
    /// MARK:: no null terminator in pattern, strlen, strdup won't work
    char * newpat = (char *)malloc(_numStrains+1);
    std::memcpy(newpat, pattern, _numStrains);
    /// MARK:: add \0 to the end
    newpat[_numStrains] = '\0';

//    if (_numDefined == _numStrains)
//        return newpat;
    /// MARK: no null terminator => strlen() = 0
    _numDefined = _numStrains;
    // remove all '?'
    unsigned int i = 0;
    while (i < _numDefined) {
        char hap = newpat[i];
        if (hap == '?') {
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
    //_pattern = newpat;
    //free(newpat);
    return newpat;
}

void MANOVA::setEigen(Dynum<std::string>& haploStrainAbbr)
{
    assert(_CorMat != nullptr);
    _useEigen = true;
    _CorMat->eigen(haploStrainAbbr);
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
    // NO null terminator in pattern, strlen, strdup won't work
    _numStrains = haploStrainAbbr.size();
    _haploStrainsAbbrevs = std::make_shared<Dynum<std::string>>(haploStrainAbbr);
    _numHaplo = this->numHaplotypes(pattern);
    // advoid memory leak
    if (_Mat != nullptr)
    {
        gsl_matrix_free(_Mat);
        _MatRowNames.clear();
        _Mat = nullptr;
        //_numDefined = 0;
    }

    if (_pattern != nullptr){
        free(_pattern);
        _pattern = nullptr;
    }
    // _pattern = pattern;
    // define _numDefined, and get new pattern
    _pattern = this->removeQMark(pattern);
    
    // Stop run if
    if (_numDefined <= _numHaplo || _numDefined < _numStrains/2 )
    {
        // std::cerr<<" numHaplo: "<<_numHaplo
        //          <<", numDefined: "<<_numDefined
        //          <<", numStrains "<<_numStrains
        //          <<", pattern: ";
        // for (unsigned int i=0; i < _numStrains; ++i) {
        //     char hap = pattern[i];
        //     if (hap != '?')
        //         std::cerr << (char) (hap + '0'); // ASCII -> char
        //     else
        //         std::cerr << hap;
        // }
        // std::cerr<<std::endl;
        return false;
    }

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
    //residuals have less rank, MANOVA cannot be performed"
    /// FIXME: _numDefined - _numHaplo == 0 
    unsigned _minL = std::min(this->L, _numDefined - _numHaplo);
    // re-assign
    //_Mat = gsl_matrix_alloc(_SubMat->size1, _SubMat->size2);
    _Mat = gsl_matrix_alloc(_numDefined, _minL);
    // skip '?' strains
    unsigned int rindex = 0, idx;
    std::string strain_abbr;
    for (size_t str1 = 0; str1 < _numStrains; str1++)
    {
        char hap = pattern[str1]; // 0,1,2,3,4, ?
        if ('?' != hap) {
            strain_abbr = _haploStrainsAbbrevs->eltOf(str1);
            //strains abbrevs exclude '?' ones
            _MatRowNames.push_back(strain_abbr);
            idx = _CorMat->eigenames.indexOf(strain_abbr);
            for (unsigned int cindex=0; cindex < _minL; ++cindex)
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
    // FIXME: need expertise to get PValue penalized.
    if (_numDefined != _numStrains) {
        FStat = INFINITY;
        PValue = 1.0;
        return;
    }
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
    for (unsigned int i = 0; i< _numHaplo; i++)
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
                 <<" pattern without qmark: "
                 <<*this<<std::endl;
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
    if (gsl_isnan(FStat))
        PValue = 1.0;
    else
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
    for(unsigned int i = 0; i < M->size1; ++i)
    {
        for (unsigned int j = 0; j< M->size2; ++j)
        {
            double x = gsl_matrix_get(M, i, j);
            gsl_matrix_set(M, i,j,std::sqrt(x));
        }
    }
}

double MANOVA::trace(const gsl_matrix *M) {
    double vv = 0;
    for (unsigned int i = 0; i< M->size1;++i)
        vv += gsl_matrix_get(M,i,i);

    return vv;
}

void MANOVA::groupMean(const gsl_matrix *M, char *pattern, gsl_matrix* _haploMean, gsl_matrix* _haploSize) {

    gsl_matrix_set_zero(_haploMean);
    gsl_matrix_set_zero(_haploSize);

    // FIXME: assert(M->size1 == strlen(_pattern) == _numDefined);
    for (unsigned int j = 0; j < M->size2; j++)
    {
        for (unsigned int i = 0; i < M->size1; i++)
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
    for (unsigned int j=0; j < M->size2; j++)
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
    for (unsigned int j=0; j < M->size2; j++) {
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

    for (unsigned int i=0; i < rows; i++){
      for(unsigned int j=0; j < cols; j++){
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
    for (unsigned int i=0; i < aov._numStrains; ++i)
    {
        char hap = aov._pattern[i];
        if (hap != '?')
            os << (char)(hap + '0'); // ASCII -> char
        else
            os << hap;
    }
    return os;
}
