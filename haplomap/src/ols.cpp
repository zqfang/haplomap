#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "ols.h"
using namespace std;

// https://github.com/cran/mvabund/blob/master/src/anova.cpp


// matrix inverse
void gsl_matrix_inv(gsl_matrix *a)
{
    size_t n=a->size1;
    //size_t m=a->size2;

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

int OLS() {
    // structure used to store Data, we want to predict
// for each person (name), the weight in function of
// gender, age, height:
// weight = c0 + c1 * gender + c2 * age + c3 * height
// we have 5 persons so N=5
// we have 4 coefficients to find: c0, c1, c2, c3, P=4
// we should obtain: -337.945  54.6301  2.96528  5.37169
    Data data[] = {
            {"Joe",	 1,	 25,	 72,	 178},
            {"Jill", 0,	 32,	 68,	 122},
            {"Jack", 1,	 27,	 69,	 167},
            {"John", 1,	 45,	 67,	 210},
            {"Jane", 0,	 38,	 62,	 108}
    };
	int N = sizeof(data) / sizeof(Data);
	int P = 4; // cste, gender, age, weight
	
	cout << "N = " << N << endl;
	cout << "P = " << P << endl;
	

	gsl_vector *y; // observed data (height)
	gsl_matrix *X; // data used to predict : cste, gender, age, weight
	gsl_vector *c; // the coefficients c0, c1, c2, c3
	gsl_matrix *cov;

	// allocate space for the matrices and vectors
	X = gsl_matrix_alloc(N, P); // this is an input
	y = gsl_vector_alloc(N); //this is an input

	c = gsl_vector_alloc(P); //this is an output
	cov = gsl_matrix_alloc(P, P); //this is an output

	//now put the data into the X matrix, row by row
	for (int i=0; i<N; ++i) {
		gsl_matrix_set(X, i, 0, static_cast<double>(1)); // because cste
		gsl_matrix_set(X, i, 1, static_cast<double>(data[i].gender));
		gsl_matrix_set(X, i, 2, static_cast<double>(data[i].age));
		gsl_matrix_set(X, i, 3, static_cast<double>(data[i].height)); 
	}
	
	// fill vector of observed data
	for (int i=0; i<N; ++i) {
		gsl_vector_set(y, i, data[i].weight);
	}
	
	double chisq;

	// allocate temporary work space for gsl 
	gsl_multifit_linear_workspace *work;
	work = gsl_multifit_linear_alloc(N, P);

	// now do the fit
	gsl_multifit_linear (X, y, c, cov, &chisq, work);

	cout << "coefficients:" << endl;
	for (int j = 0; j<P; j++) {
		cout << "c" << j << " = " << std::setprecision(9);
		cout << gsl_vector_get(c, j) << endl;
	}
	
	cout << endl;
	cout << "expected <=> predicted" << endl;
	for (int i=0; i<N; ++i) {
		double r = gsl_vector_get(c, 0);
		r += data[i].gender * gsl_vector_get(c, 1); 
		r += data[i].age * gsl_vector_get(c, 2);
		r += data[i].height * gsl_vector_get(c, 3);
		cout << data[i].weight << " <=> " << std::setprecision(9) << r << endl;
	}
	
	//************************************
	//very important stuff at the end here
	//************************************
	// don't forget to deallocate - there is no garbage collector in C!

	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(c);

	return 0;
}
