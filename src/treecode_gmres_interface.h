/**************************************************************************
* FILE NAME: treecode_gmres_interface.h                                   *
*                                                                         *
* PURPOSE: header for functions that interface treecode with GMRes,       *
*          namely, the matrix multiply and solve functions                *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
**************************************************************************/

#ifndef H_TREECODE_GMRES_INTERFACE_H
#define H_TREECODE_GMRES_INTERFACE_H

#include "struct_particles.h"

/* functions used by GMRes */
int matvec(double *alpha, double *tpoten_old,
           double *beta, double *tpoten, struct Particles *particles);

int psolve(double *z, double *r, struct Particles *particles);

int psolve_precond(double *z, double *r, struct Particles *particles);

#endif /* H_TREECODE_GMRES_INTERFACE_H */
