/*
 * Modified from original CLAPACK version to support MKL.
 * Original file description below.
 */

/*                                                                                           
 * This file is part of CLAPACK, the f2c'ed version of LAPACK                                
 * hosted at netlib.org.                                                                     
 *                                                                                           
 * For more information, see the LAPACK User's Guide:                                        
 *                                                                                           
  @BOOK{laug,                                                                                
        AUTHOR = {Anderson, E. and Bai, Z. and Bischof, C. and                               
                  Blackford, S. and Demmel, J. and Dongarra, J. and                          
                  Du Croz, J. and Greenbaum, A. and Hammarling, S. and                       
                  McKenney, A. and Sorensen, D.},                                            
        TITLE = {{LAPACK} Users' Guide},                                                     
        EDITION = {Third},                                                                   
        PUBLISHER = {Society for Industrial and Applied Mathematics},                        
        YEAR = {1999},                                                                       
        ADDRESS = {Philadelphia, PA},                                                        
        ISBN = {0-89871-447-8 (paperback)}  }                                                
 *                                                                                           
 */ 

/* GMRES.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
 */

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#include <cmath>
#include <cstdio>
#include <mkl_cblas.h>

#include "treecode.h"

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  GMRES solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  RESTRT  (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix H (see WORK and H).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,RESTRT+4).
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  H       (workspace) DOUBLE PRECISION array, dimension (LDH,RESTRT+2).
*          This workspace is used for constructing and storing the
*          upper Hessenberg matrix. The two extra columns are used to
*          store the Givens rotation matrices.
*
*  LDH    (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,RESTRT+1).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common
*          block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: LDH < RESTRT
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DROT, DROTG, DSCAL
*  ============================================================
*/

static int update_(long int, long int, double *, double *, long int,
                   double *, double *, double *, long int);
static  int basis_(long int, long int, double *, double *, long int, double *);

//*****************************************************************
int Treecode::gmres_(long int n, const double *b, double *x,
                     long int restrt, long int ldw, long int ldh,
                     long int& iter, double& resid)
{
    std::vector<double> work_vec(ldw * (restrt + 4));
    std::vector<double> h_vec   (ldh * (restrt + 2));
    
    double* work = work_vec.data();
    double* h    = h_vec.data();


/*     Test the input parameters. */

    if (n < 0) {
        return -1;
    } else if (ldw < MAX(1,n)) {
        return -2;
    } else if (iter <= 0) {
        return -3;
    } else if (ldh < restrt + 1) {
        return -4;
    }

    long int maxit = iter;
    double tol = resid;


/*     Store the Givens parameters in matrix H. */
/*     Set initial residual (AV is temporary workspace here). */

    cblas_dcopy(n, b, 1, &work[2 * ldw], 1);

    if (cblas_dnrm2(n, x, 1) != 0.) {
        cblas_dcopy(n, b, 1, &work[2 * ldw], 1);
        Treecode::matrix_vector(-1., x, 1., &work[2 * ldw]);
    }

    Treecode::precondition(work, &work[2 * ldw]);

    double bnrm2 = cblas_dnrm2(n, b, 1);
    if (bnrm2 == 0.) {
        bnrm2 = 1.;
    }
    
    if (cblas_dnrm2(n, work, 1) / bnrm2 < tol) {
	    return 0;
    }

    iter = 0;

    while (true) {

    /*        Construct the first column of V. */

        cblas_dcopy(n, work, 1, &work[3 * ldw], 1);
        double rnorm = cblas_dnrm2(n, &work[3 * ldw], 1);
        cblas_dscal(n, 1. / rnorm, &work[3 * ldw], 1);

    /*        Initialize S to the elementary vector E1 scaled by RNORM. */

        work[ldw] = rnorm;
        for (long int k = 1; k < n; ++k) {
            work[k + ldw] = 0.;
        }

        for (long int i = 0; i < restrt; ++i) {
            ++iter;

            Treecode::matrix_vector(1., &work[(3 + i) * ldw], 0., &work[2 * ldw]);
            
            Treecode::precondition(&work[2 * ldw], &work[2 * ldw]);

        /*           Construct I-th column of H orthnormal to the previous */
        /*           I-1 columns. */

            basis_(i+1, n, &h[i * ldh], &work[3 * ldw], ldw, &work[2 * ldw]);

        /*           Apply Givens rotations to the I-th column of H. This */
        /*           "updating" of the QR factorization effectively reduces */
        /*           the Hessenberg matrix to upper triangular form during */
        /*           the RESTRT iterations. */

            for (long int k = 0; k < i; ++k) {
                cblas_drot(1, &h[k + i * ldh], ldh, &h[k + 1 + i * ldh],
                           ldh, h[k + restrt * ldh], h[k + (restrt + 1) * ldh]);
            }

        /*           Construct the I-th rotation matrix, and apply it to H so that */
        /*           H(I+1,I) = 0. */

            double aa = h[i * (ldh + 1)];
            double bb = h[i * (ldh + 1) + 1];
            
            cblas_drotg(&aa, &bb, &h[i + (restrt)     * ldh],
                                  &h[i + (restrt + 1) * ldh]);
                                  
            cblas_drot(1,   &h[i * (ldh + 1)],
                       ldh, &h[i * (ldh + 1) + 1],
                       ldh,  h[i + (restrt)     * ldh],
                             h[i + (restrt + 1) * ldh]);

        /*           Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This */
        /*           gives an approximation of the residual norm. If less than */
        /*           tolerance, update the approximation vector X and quit. */

            cblas_drot(1,   &work[i + ldw],
                       ldw, &work[i + ldw + 1],
                       ldw, h[i + (restrt)     * ldh],
                            h[i + (restrt + 1) * ldh]);
                            
            resid = std::fabs(work[i + 1 + ldw]) / bnrm2;

            std::printf("iteration no. = %ld, error = %e\n", iter, resid);

            if (resid <= tol) {
                update_(i+1, n, x, h, ldh, &work[2 * ldw], &work[ldw], &work[3 * ldw], ldw);
                return 0;
            }
        }

    /*        Compute current solution vector X. */

        update_(restrt, n, x, h, ldh, &work[2 * ldw], &
                work[ldw], &work[3 * ldw], ldw);

    /*        Compute residual vector R, find norm, then check for tolerance. */

        cblas_dcopy(n, b, 1, &work[2 * ldw], 1);
        
        Treecode::matrix_vector(-1., x, 1., &work[2 * ldw]);
        Treecode::precondition(work, &work[2 * ldw]);
        
        work[restrt + ldw] = cblas_dnrm2(n, work, 1);
        resid = work[restrt + ldw] / bnrm2;
        
        if (resid <= tol) {
            return 0;
        }
        
        if (iter == maxit) {
            return 1;
        }

/*        Restart. */
    }
}


/*     =============================================================== */
static int update_(long int i, long int n, double *x, double *h, long int ldh,
                   double *y, double *s, double *v, long int ldv)
{
/*     This routine updates the GMRES iterated solution approximation. */

/*     Solve H*Y = S for upper triangualar H. */

    cblas_dcopy(i, s, 1, y, 1);
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
                i, h, ldh, y, 1);

/*     Compute current solution vector X = X + V*Y. */

    cblas_dgemv(CblasColMajor, CblasNoTrans, n, i, 1.,
                v, ldv, y, 1, 1., x, 1);

    return 0;
}


/*     ========================================================= */
static int basis_(long int i, long int n, double *h, double *v, long int ldv, double *w)
{
/*     Construct the I-th column of the upper Hessenberg matrix H */
/*     using the Gram-Schmidt process on V and W. */

    for (long int k = 0; k < i; ++k) {
        h[k] = cblas_ddot(n, w, 1, &v[k * ldv], 1);
        cblas_daxpy(n, -h[k], &v[k * ldv], 1, w, 1);
    }
    
    h[i] = cblas_dnrm2(n, w, 1);
    cblas_dcopy(n, w, 1, &v[i * ldv], 1);
    cblas_dscal(n, 1. / h[i], &v[i * ldv], 1);

    return 0;
}
