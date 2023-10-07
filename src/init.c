#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(em_diag)(double *longx, int *xdims, int *ndim, double *y, int *n, 
int *nclass, int *subnum, double *longmu, double *diagsig_inv_sq, double *pis, 
double *longgamma, int *glen, int *maxiter, double *loglike);
extern void F77_NAME(pred_diag)(double *longx, int *xdims, int *ndim, double *prob, 
int *n, int *nclass, int *subnum, double *longmu, double *diagsig_inv_sq, 
double *pis, double *phi);
extern void F77_NAME(em)(double *longx, int *xdims, int *ndim, double *y, int *n, 
int *nclass, int *subnum, double *longmu, double *longsig_inv_sq, double *pis, 
double *longgamma, int *glen, int *maxiter, double *loglike);
extern void F77_NAME(pred_f)(double *longx, int *xdims, int *ndim, double *prob, 
int *n, int *nclass, int *subnum, double *longmu, double *longsig_inv_sq, 
double *pis, double *phi); 
extern void F77_NAME(tucker_prod_f)(double *tnew, int *newdims, double *tnsr, 
int *tnsrdims, int *ndim, double *longvec, int *veclen, int *skip);

static const R_FortranMethodDef FortranEntries[] = {
    {"em_diag", (DL_FUNC) &F77_NAME(em_diag), 14},
    {"pred_diag", (DL_FUNC) &F77_NAME(pred_diag), 11},
    {"em", (DL_FUNC) &F77_NAME(em), 14},
    {"pred_f", (DL_FUNC) &F77_NAME(pred_f), 11},
    {"tucker_prod_f", (DL_FUNC) &F77_NAME(tucker_prod_f), 8},
    {NULL, NULL, 0}
};

void R_init_tmda(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
