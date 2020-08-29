#ifndef H_GMRES_H
#define H_GMRES_H

extern "C" {
int gmres_(long int n, double *b, double *x, long int *restrt, double *work, 
           long int ldw, double *h, long int ldh, long int *iter, double *resid, 
           int (*matvec) (double *, double *, double *, double *), 
           int (*psolve) (double *, double *), 
           long int *info);
}

#endif
