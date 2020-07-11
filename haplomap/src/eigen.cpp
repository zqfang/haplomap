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

EigenMat::EigenMat(const char* MatrixFile, const char* RowNamesFile)
{
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
    // read matrix file
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

    //    std::ifstream input(filename);
    //    string line;
    //    if (input.is_open()) {
    //        while (getline(input, line))
    //            std::cout << line <<'\n';
    //    }
    //    input.close();
}
EigenMat::EigenMat(const std::vector<std::vector<double>> &Mat)
{
    size1 = Mat.size();
    size2 = Mat[0].size();
    data = gsl_matrix_alloc(size1, size2);
    for (int i = 0; i < data->size1; ++i){
        for (int j=0; j < data->size2; ++j){
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
        gsl_vector_free(eigenvalues);
        gsl_matrix_free(eigenvectors);
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

void EigenMat::calcVariance() {
    if (data == nullptr) {
        std::cerr<<"Error! data not found!";
        return;
    }

    if (eigenvectors == nullptr) {
        this->eigen();
    }
    if (variances != nullptr)
        return;

    variances = gsl_vector_alloc(size1);
    gsl_vector_memcpy(variances, eigenvalues);

    double _sumEigen = 0;
    for (int i = 0; i < size1; ++i)
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
    gsl_vector* mean = gsl_vector_alloc(rows);

    for(i = 0; i < rows; i++) {
        gsl_vector_set(mean, i, gsl_stats_mean(data->data + i * cols, 1, cols));
    }

    // Get mean-substracted data into matrix mean_substracted_data.
    gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
    gsl_matrix_memcpy(mean_substracted_data, data);
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
    gsl_vector* _eigenvalues = gsl_vector_alloc(rows);
    gsl_matrix* _eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(covariance_matrix, _eigenvalues, _eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(covariance_matrix);

    // Sort the eigenvectors
    gsl_eigen_symmv_sort(_eigenvalues, _eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
    gsl_vector_free(_eigenvalues);

    // Project the original dataset
    gsl_matrix* result = gsl_matrix_alloc(L, cols);
    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(_eigenvectors, 0, 0, rows, L);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, data, 0.0, result);
    gsl_matrix_free(_eigenvectors);
    // Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L
    return result;
}