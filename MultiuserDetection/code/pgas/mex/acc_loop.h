#include <math.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
//#include "gsl_headers.h"
//#include "matrix.h"


#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
//#include <gsl/gsl_math.h>
//#include "gsl/gsl_cdf.h"
//#include "gsl/gsl_randist.h"

#include "mex.h"

//Methods included in this header file
inline double get_Y(double *Y, int n, int t, int Nr, int T);
inline double get_H(double *H, int nr, int nt, int ll, int Nr, int Nt, int L);
inline int getint_2D(double *x, int n1, int n2, int N1, int N2);
inline int getint_3D(double *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline void setdouble_2D(double val, double *x, int n1, int n2, int N1, int N2);
inline double getdouble_2D(double *x, int n1, int n2, int N1, int N2);
