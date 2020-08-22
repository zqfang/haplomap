//
// Created by Zhuoqing Fang on 7/8/20.
//

#include <assert.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include "gsl/gsl_matrix.h"
#include "ColumnReader.h"
#include "eigen.h"

namespace HBCGM {
    // debugging
    void print_vec(std::vector<double> &M)
    {
        std::cout<<"trace vector: "<<std::endl;
        for (auto &v: M) {
            std::cout<<v<<" ";
        }
        std::cout<<std::endl;
    }
    void print_gvec(gsl_vector* M)
    {
        std::cout<<"trace gsl vector: "<<std::endl;
        for (unsigned i = 0; i < M->size; ++i)
        {
            std::cout<<gsl_vector_get(M, i)<<" ";
        }
        std::cout<<std::endl;
    }

    void print_gmat(gsl_matrix* M){
        std::cout << "trace gsl matrix: "<<std::endl;
        for (unsigned i = 0; i < M->size1; ++i){
            for (unsigned j=0; j < M->size2; ++j) {
                std::cout<< gsl_matrix_get(M, i, j) << " ";
            }
            std::cout<<std::endl;
        }
    }
}


EigenMat::EigenMat(const char* MatrixFile, const char* RowNamesFile){
    // read row names file
    ColumnReader rdr(RowNamesFile, (char *)"\t");
    int numtoks;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        std::string strain_abbrev = rdr.getToken(1);
        int strIdx = rownames.addElementIfNew(strain_abbrev);
        if (strIdx < 0)
        {
            std::cout << "Undefined strain abbrev: " << strain_abbrev << std::endl;
        }
    }
    ColumnReader rmat(MatrixFile, (char *)"\t");
    data = gsl_matrix_alloc(rownames.size(), rownames.size());
    int i = 0;
    while ((numtoks = rmat.getLine()) >= 0) {
        for (int t = 0; t < numtoks; t++){
            double d = std::stod(rmat.getToken(t));
            gsl_matrix_set(data,i,t,d);
        }
        i++;
    }
    size1 = data->size1;
    size2 = data->size2;
}

EigenMat::EigenMat(const char* filename, bool header, const std::string delimimiter )
{
    int numtoks;
    std::vector<std::vector<double>> temp_mat;
    // read matrix file
    ColumnReader rmat(filename, (char *)delimimiter.c_str());
    // get first line
    numtoks = rmat.getLine();
    if (numtoks < 0)
        std::invalid_argument("Empty input file");

    std::vector<double> temp;
    // read header
    if (header) {
        for (int t = 0; t < numtoks; t++)
            rownames.addElementIfNew(rmat.getToken(t));
    } else {
        for (int t = 0; t < numtoks; t++)
            temp.push_back(std::stod(rmat.getToken(t)));
        temp_mat.push_back(temp);
        temp.clear();
    }
    // read the rest of file
    while ((numtoks = rmat.getLine()) >= 0) {
        temp.clear();
        for (int t = 0; t < numtoks; t++) {
            temp.push_back(std::stod(rmat.getToken(t)));
        }
        temp_mat.push_back(temp);
    }
    size1 = temp_mat.size();
    size2 = temp_mat[0].size();
    this->data = gsl_matrix_alloc(size1, size2);
    for (unsigned i = 0; i < size1; ++i) {
        for (unsigned j = 0; j < size2; ++j) {
            gsl_matrix_set(this->data, i, j, temp_mat[i][j]);
            //std::cout << temp_mat[i][j] << " ";
        }
        //std::cout<<std::endl;
    }
}

EigenMat::EigenMat(const std::vector<std::vector<double>> &Mat)
{
    size1 = Mat.size();
    size2 = Mat[0].size();
    data = gsl_matrix_alloc(size1, size2);
    for (size_t i = 0; i < data->size1; ++i){
        for (size_t j=0; j < data->size2; ++j){
            gsl_matrix_set(data,i,j,Mat[i][j]);
        }
    }
}

EigenMat::EigenMat(const gsl_matrix * Mat)
{
    data = gsl_matrix_alloc(Mat->size1, Mat->size2);
    gsl_matrix_memcpy(data, Mat);
    size1 = Mat->size1;
    size2 = Mat->size2;
}

EigenMat::~EigenMat()
{
    gsl_matrix_free(data);
    if (eigenvectors != nullptr)
    {
        gsl_vector_free(eigenvalues);
        gsl_matrix_free(eigenvectors);
    }
    if (variances != nullptr)
        gsl_vector_free(variances);
}

void EigenMat::eigen(){
    if (data == nullptr) {
        std::cerr<<"Error! data not found!";
        return;
    }
    // already calculated
    if (eigenvectors != nullptr)
        return;

    // calc
    unsigned int rows = data->size1;
    // unsigned int cols = data->size2;

    // Get eigenvectors, sort by eigenvalue.
    eigenvalues = gsl_vector_alloc(rows);
    eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(data, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    // Sort the eigenvectors
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
}

void EigenMat::eigen(Dynum<std::string> &strainAbbrev) {
    if (data == nullptr) {
        std::cerr<<"Error! data not found!";
        return;
    }

    eigenames = strainAbbrev;
    // already calculated
    if (eigenvectors != nullptr)
        return;

     std::string strain_abbr;
     int newind;
     gsl_matrix *subdata = gsl_matrix_alloc(strainAbbrev.size(), strainAbbrev.size());
     for (int str1=0; str1 < strainAbbrev.size(); str1++){
         strain_abbr = strainAbbrev.eltOf(str1);
         newind = rownames.indexOf(strain_abbr);
         for (int cind = 0; cind < strainAbbrev.size(); cind ++)
             gsl_matrix_set(subdata, str1, cind, gsl_matrix_get(data, newind, cind));
     }
    // calc
    unsigned int rows = subdata->size1;
    // unsigned int cols = data->size2;
    // Get eigenvectors, sort by eigenvalue.
    eigenvalues = gsl_vector_alloc(rows);
    eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(subdata, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    // Sort the eigenvectors
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
    gsl_matrix_free(subdata);
}

void EigenMat::calcVariance() {
    if (data == nullptr) {
        std::cerr<<"Error! data not found!";
        return;
    }

    if (eigenvectors == nullptr)
    {
        this->eigen();
    }
    if (variances != nullptr)
        return;

    variances = gsl_vector_alloc(eigenvalues->size);
    gsl_vector_memcpy(variances, eigenvalues);

    double _sumEigen = 0;
    for (unsigned i = 0; i < eigenvalues->size; ++i)
        _sumEigen += gsl_vector_get(variances, i);
    gsl_vector_scale(variances, 100.0/_sumEigen);
}

/*
@param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
@param L - dimension reduction
*/
gsl_matrix* EigenMat::pca(unsigned int L)
{
    assert(data != NULL);
    assert(L > 0 && L < data->size2);

    unsigned int i;
    unsigned int rows = data->size1;
    unsigned int cols = data->size2;
    //HBCGM::print_gmat(data);
    gsl_vector* mean = gsl_vector_alloc(rows);
    // FIXME: why could not allocate memory here
    // get row means
    for(i = 0; i < rows; i++) {
        //get current column as vector
        gsl_vector_view myRow = gsl_matrix_row(this->data, i);
        // calc mean
        double avg = gsl_stats_mean(myRow.vector.data, myRow.vector.stride, myRow.vector.size);
        // assign
        gsl_vector_set(mean, i, avg);

        //gsl_vector_set(mean, i, gsl_stats_mean(data->data + i * cols, 1, cols));
    }

    // Get mean-substracted data into matrix mean_substracted_data.
    gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
    gsl_matrix_memcpy(mean_substracted_data, this->data);
    for(i = 0; i < cols; i++) {
        gsl_vector_view mean_substracted_point_view = gsl_matrix_column(mean_substracted_data, i);
        gsl_vector_sub(&mean_substracted_point_view.vector, mean);
    }
    gsl_vector_free(mean);

    // Compute Covariance matrix
    gsl_matrix* covariance_matrix = gsl_matrix_alloc(rows, rows);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0 / (double)(cols - 1),
            mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
    gsl_matrix_free(mean_substracted_data);

    // Get eigenvectors, sort by eigenvalue.
    eigenvalues = gsl_vector_alloc(rows);
    eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(covariance_matrix);

    // Sort the eigenvectors
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

    // Project the original dataset
    gsl_matrix* result = gsl_matrix_alloc(L, cols);
    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, data, 0.0, result);
    // Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L
    return result;
}
