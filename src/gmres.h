#ifndef H_GMRES_H
#define H_GMRES_H

#include "struct_particles.h"

int gmres_(long int n, double *b, double *x, long int *restrt, double *work, 
           long int ldw, double *h, long int ldh, long int *iter, double *resid, 
           int (*matvec) (double *, double *, double *, double *, struct Particles *),
           int (*psolve) (double *, double *, struct Particles *), 
           long int *info, struct Particles *particles);

#endif
