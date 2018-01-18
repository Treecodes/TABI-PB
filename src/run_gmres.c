/**************************************************************************
* FILE NAME: run_gmres.c                                                  *
*                                                                         *
* PURPOSE: Wraps GMRes routine called by tabipb.c, links to psolve and    *
*          matvec routines located in treecode.c                          *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 01/11/2018  Leighton Wilson   Created, moved from tabipb routine        *
*                                                                         *
**************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "run_gmres.h"
#include "treecode_gmres_interface.h"

#include "array.h"
#include "particle_struct.h"

/********************************************************/
int RunGMRES(int nface, double *source_term, double *xvct)
{
    int i;
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
    
    for (i = 0; i < N; i++) xvct[i] = 0.0;
    
    gmres_(&N, source_term, xvct, &RESTRT, work, &ldw, 
           h, &ldh, &iter, &resid, &matvec, psolve, &info);

    free_vector(work);
    free_vector(h);

    return 0;
}
/********************************************************/
