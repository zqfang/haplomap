#include <assert.h>
#include <gsl/gsl_blas.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_complex_math.h> 
#include <gsl/gsl_linalg.h>

gsl_matrix* pca(const gsl_matrix* data, unsigned int L)
{
    /*
    @param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
    @param L - dimension reduction
    */
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
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
    gsl_matrix_free(mean_substracted_data);

    // Get eigenvectors, sort by eigenvalue.
    gsl_vector* eigenvalues = gsl_vector_alloc(rows);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(covariance_matrix);

    // Sort the eigenvectors
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
    gsl_vector_free(eigenvalues);

    // Project the original dataset
    gsl_matrix* result = gsl_matrix_alloc(L, cols);
    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, data, 0.0, result);
    gsl_matrix_free(eigenvectors);

    // Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L
    return result;
}

// matrix multiplication GSL_BLAS
// https://www.gnu.org/software/gsl/doc/html/blas.html
// s,d,c, z = double precision ...
// ge = general matrices 
// mm = matrix-matrix multiplication 
// gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &B.matrix, 0.0, &C.matrix); 
// example
gsl_matrix_complex *multiply(gsl_matrix_complex *A, gsl_matrix_complex *B) 
{ 
    gsl_matrix_complex *result = gsl_matrix_complex_alloc(A->size1, B->size2); 

    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 
        GSL_COMPLEX_ONE, A, B, 
        GSL_COMPLEX_ZERO, result); 

    return result; 
} 

// or
// matrix multiplication 
// c=a*b, a: M*N,b: N*L, then c: M*L.
void gsl_matrix_mul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c)
{
	for (size_t i=0;i<a->size1;i++)
	{
		for (size_t j=0;j<b->size2;j++)
		{
			double sum=0.0;
			for (size_t k=0;k<b->size1;k++)
			{
				sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
			}
			gsl_matrix_set(c,i,j,sum);
		}
	}
}

// matrix inverse
void gsl_matrix_inv(gsl_matrix *a)
{
	size_t n=a->size1;
	size_t m=a->size2;
 
	gsl_matrix *temp1=gsl_matrix_calloc(n,n);
	gsl_matrix_memcpy(temp1,a);
 
	gsl_permutation *p=gsl_permutation_calloc(n);
	int sign=0;
	gsl_linalg_LU_decomp(temp1,p,&sign);
	gsl_matrix *inverse=gsl_matrix_calloc(n,n);
 
	gsl_linalg_LU_invert(temp1,p,inverse);
	gsl_matrix_memcpy(a,inverse);
 
	gsl_permutation_free(p);
	gsl_matrix_free(temp1);
	gsl_matrix_free(inverse);
 
}