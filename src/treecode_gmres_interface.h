/*
 * C header file for treecode GMRes interface routines of tabipb
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/12/2018
 */

#ifndef H_TREECODE_GMRES_INTERFACE_H
#define H_TREECODE_GMRES_INTERFACE_H


/* functions used by GMRes */
int matvec(double *alpha, double *tpoten_old,
           double *beta, double *tpoten);

int psolve(double *z, double *r);

#endif /* H_TREECODE_GMRES_INTERFACE_H */
