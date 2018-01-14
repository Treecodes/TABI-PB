/*
 * C routine for setting up and calling GMRes
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/11/2018
 */

#include <math.h>
#include <stdlib.h>

#include "run_gmres.h"
#include "treecode_gmres_interface.h"

#include "array.h"
#include "particle_struct.h"


int RunGMRES(int nface, double *source_term, double *xvct)
{
    static long int info;
    long int RESTRT, ldw, ldh, iter, N;
    double resid;
    double *work, *h;

    extern int gmres_(long int *n, double *b, double *x, long int *restrt,
                      double *work, long int *ldw, double *h, long int *ldh,
                      long int *iter, double *resid, int (*matvec)(),
                      int (*psolve)(), long int *info);

    /* parameters for GMRES */
    RESTRT = 10;
    N = 2 * nface;
    ldw = N;
    ldh = RESTRT + 1;
    iter = 100;
    resid = 1e-4;

    make_vector(work, ldw * (RESTRT + 4));
    make_vector(h, ldh * (RESTRT + 2));

    gmres_(&N, source_term, xvct, &RESTRT, work, &ldw, 
           h, &ldh, &iter, &resid, &matvec, psolve, &info);

    free_vector(work);
    free_vector(h);


    return 0;
}
