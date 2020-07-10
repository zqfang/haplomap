#include <assert.h>
#include <gsl/gsl_blas.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_complex_math.h> 
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>


gsl_matrix* pca(const gsl_matrix* data, unsigned int L);

// matrix multiplication GSL_BLAS
// https://www.gnu.org/software/gsl/doc/html/blas.html
// s,d,c, z = double precision ...
// ge = general matrices 
// mm = matrix-matrix multiplication 
// gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &B.matrix, 0.0, &C.matrix); 
// example
gsl_matrix_complex *multiply(gsl_matrix_complex *A, gsl_matrix_complex *B);

// or
// matrix multiplication 
// c=a*b, a: M*N,b: N*L, then c: M*L.
void gsl_matrix_matmul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c);

// matrix inverse
void gsl_matrix_inv(gsl_matrix *a);


//// get eigen values
//gsl_matrix getEigenValues(gsl_matrix & M) {
//    int k = M.ncol();
//
//    RcppGSL::Vector ev(k);  	// instead of gsl_vector_alloc(k);
//    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(k);
//    gsl_eigen_symm(M, ev, w);
//    gsl_eigen_symm_free (w);
//
//    return ev;				// return results vector
//}

typedef struct {
    std::string name;
    int gender; // 1 male, 0 female
    int age;    // years
    int height; // inches
    int weight; // pounds
} Data;