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
**************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "run_gmres.h"
#include "treecode_gmres_interface.h"
#include "gmres.h"

#include "array.h"
#include "struct_particles.h"

/********************************************************/
int RunGMRES(int nface, double *source_term, int precond,
             double *xvct, long int *iter, struct Particles *particles)
{
    static long int info;

    /* parameters for GMRES */
    long int RESTRT = 10;
    long int N = 2 * nface;
    long int ldw = N;
    long int ldh = RESTRT + 1;
    double resid = 1e-4;
        
    *iter = 100;

    double *work = (double *)malloc(ldw * (RESTRT + 4) * sizeof(double));
    double *h    = (double *)malloc(ldh * (RESTRT + 2) * sizeof(double));
    
    for (long int i = 0; i < N; i++) xvct[i] = 0.0;
    
    if (precond == 0) {
        gmres_(N, source_term, xvct, &RESTRT, work, ldw,
               h, ldh, iter, &resid, &matvec, psolve, &info, particles);
    } else if (precond == 1) {
        gmres_(N, source_term, xvct, &RESTRT, work, ldw,
               h, ldh, iter, &resid, &matvec, psolve_precond, &info, particles);
    }

    free(work);
    free(h);

    return 0;
}
/********************************************************/
