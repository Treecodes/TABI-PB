#ifndef H_GMRES_BLAS_H
#define H_GMRES_BLAS_H

#include "f2c.h"

extern "C" {
int drot_(integer n, doublereal *dx, integer incx,
    doublereal *dy, integer incy, doublereal c, doublereal s);
 
int drotg_(doublereal *da, doublereal *db, doublereal *c,
    doublereal *s);

int dtrsv_(char *uplo, char *trans, char *diag, integer n,
    doublereal *a, integer lda, doublereal *x, integer incx);
 
int dgemv_(char *trans, integer m, integer n, doublereal alpha,
    doublereal *a, integer lda, doublereal *x, integer incx,
    doublereal beta, doublereal *y, integer incy);
}
#endif
