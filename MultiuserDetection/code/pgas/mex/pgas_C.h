#include <math.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
//#include "gsl_headers.h"
#include <matrix.h>


#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
//#include <gsl/gsl_math.h>
//#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

#include "mex.h"


#define myINFINITY 1.0e100

//Methods included in this header file
inline double getdouble_2D(double *x, int n1, int n2, int N1, int N2);
inline double getdouble_3D(double *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline void setdouble_2D(double val, double *x, int n1, int n2, int N1, int N2);
inline void setdouble_3D(double val, double *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline int getint_2D(double *x, int n1, int n2, int N1, int N2);
inline int getint_2D(int *x, int n1, int n2, int N1, int N2);
inline int getint_forC_2D(int *x, int n1, int n2, int N1, int N2, double *head, int offState);
inline int getint_3D(double *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline int getint_3D(int *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline int getint_Xt(int16_T *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline int getint_forC_Xt(int16_T *x, int n1, int n2, int n3, int N1, int N2, int N3, double *head, int offState);
inline void setint_2D(int val, int *x, int n1, int n2, int N1, int N2);
inline void setint_3D(int val, int *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline void setint_Xt(int val, int16_T *x, int n1, int n2, int n3, int N1, int N2, int N3);
inline bool getbool_2D(bool *x, int n1, int n2, int N1, int N2);
inline void setbool_2D(int val, bool *x, int n1, int n2, int N1, int N2);

// Misc functions
inline double my_exp(double val);
inline double my_log(double val);

// Auxiliary functions
double *computeFactor1(int Nt, int Nr, int N, int t, int T, int L, int Q, double *Cr, double *Ci, double *Yr, double *Yi, int *a_ind, int *xc, int16_T *Xt, double *Hr, double *Hi, bool flagComplex, mxComplexity myComplex, double sy2, double *head, int offState);
double *computeFactor2(int Nt, int N, int t, int T, int Q, double logQ_1, int *xc, int16_T *Xt, double *logam, double *logbm, double *log1_am, double *log1_bm, int lHead, int onOffModel, int offState);
