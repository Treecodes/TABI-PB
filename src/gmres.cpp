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

#define max(a,b) ((a) >= (b) ? (a) : (b))

#include <stdio.h>
#include <math.h>
#include <mkl_cblas.h>

#include "struct_particles.h"

#ifdef MPI_ENABLED
    #include <mpi.h>
#endif

/* Table of constant values */

static long int c__1 = 1;
static double c_b7 = -1.;
static double c_b8 = 1.;
static double c_b20 = 0.;

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


//*****************************************************************
int gmres_(long int n, double *b, double *x, long int *restrt, double *work, 
           long int ldw, double *h, long int ldh, long int *iter, double *resid, 
           int (*matvec) (double *, double *, double *, double *, struct Particles *),
           int (*psolve) (double *, double *, struct Particles *), long int *info,
           struct Particles *particles)
{
    /* System generated locals */
    long int work_offset, h_offset, i__1;
    double d__1;

    /* Local variables */
    static double bnrm2, rnorm, aa, bb, tol;
    static long int i, k, r, s, v, w, y, maxit, cs, av, sn;
    extern /* Subroutine */ int update_(long int *, long int, double *, double *, long int,
                                        double *, double *, double *, long int);
    extern /* Subroutine */ int basis_(long int *, long int, double *, double *, long int, double *);
    
    int rank = 0;
    
#ifdef MPI_ENABLED
    int ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    /* Parameter adjustments */

    h_offset = ldh + 1;
    h -= h_offset;
    work_offset = ldw + 1;
    work -= work_offset;
    --x;
    --b;

    /* Executable Statements */
    *info = 0;

/*     Test the input parameters. */

    if (n < 0) {
	*info = -1;
    } else if (ldw < max(1,n)) {
	*info = -2;
    } else if (*iter <= 0) {
	*info = -3;
    } else if (ldh < *restrt + 1) {
	*info = -4;
    }
    if (*info != 0) {
	return 0;
    }

    maxit = *iter;
    tol = *resid;

/*     Alias workspace columns. */

    r = 1;
    s = r + 1;
    w = s + 1;
    y = w;
    av = y;
    v = av + 1;

/*     Store the Givens parameters in matrix H. */

    cs = *restrt + 1;
    sn = cs + 1;

/*     Set initial residual (AV is temporary workspace here). */

    cblas_dcopy(n, &b[1], c__1, &work[av * ldw + 1], c__1);

    if (cblas_dnrm2(n, &x[1], c__1) != 0.) {
/*        AV is temporary workspace here. */

		cblas_dcopy(n, &b[1], c__1, &work[av * ldw + 1], c__1);
		(*matvec)(&c_b7, &x[1], &c_b8, &work[av * ldw + 1], particles);
	}

    (*psolve)(&work[r * ldw + 1], &work[av * ldw + 1], particles);

    bnrm2 = cblas_dnrm2(n, &b[1], c__1);
    if (bnrm2 == 0.) {
	bnrm2 = 1.;
    }
    if (cblas_dnrm2(n, &work[r * ldw + 1], c__1) / bnrm2 < tol) {
	goto L70;
    }

    *iter = 0;

L10:

    i = 0;

/*        Construct the first column of V. */

    cblas_dcopy(n, &work[r * ldw + 1], c__1, &work[v * ldw + 1], c__1);
    rnorm = cblas_dnrm2(n, &work[v * ldw + 1], c__1);
    d__1 = 1. / rnorm;
    cblas_dscal(n, d__1, &work[v * ldw + 1], c__1);

/*        Initialize S to the elementary vector E1 scaled by RNORM. */

    work[s * ldw + 1] = rnorm;
    i__1 = n;
    for (k = 2; k <= i__1; ++k) {
	work[k + s * ldw] = 0.;
/* L20: */
    }

L30:

    ++i;
    ++(*iter);

    (*matvec)(&c_b8, &work[(v + i - 1) * ldw + 1], &c_b20, &work[av *
	    ldw + 1], particles);
    
	(*psolve)(&work[w * ldw + 1], &work[av * ldw + 1], particles);

/*           Construct I-th column of H orthnormal to the previous */
/*           I-1 columns. */

    basis_(&i, n, &h[i * ldh + 1], &work[v * ldw + 1], ldw, &work[w *
	     ldw + 1]);

/*           Apply Givens rotations to the I-th column of H. This */
/*           "updating" of the QR factorization effectively reduces */
/*           the Hessenberg matrix to upper triangular form during */
/*           the RESTRT iterations. */

    i__1 = i - 1;
    for (k = 1; k <= i__1; ++k) {
	cblas_drot(c__1, &h[k + i * ldh], ldh, &h[k + 1 + i * ldh], ldh, h[
		k + cs * ldh], h[k + sn * ldh]);
/* L40: */
    }

/*           Construct the I-th rotation matrix, and apply it to H so that
 */
/*           H(I+1,I) = 0. */

    aa = h[i + i * ldh];
    bb = h[i + 1 + i * ldh];
    cblas_drotg(&aa, &bb, &h[i + cs * ldh], &h[i + sn * ldh]);
    cblas_drot(c__1, &h[i + i * ldh], ldh, &h[i + 1 + i * ldh], ldh, h[i + 
	    cs * ldh], h[i + sn * ldh]);

/*           Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This */
/*           gives an approximation of the residual norm. If less than */
/*           tolerance, update the approximation vector X and quit. */

    cblas_drot(c__1, &work[i + s * ldw], ldw, &work[i + 1 + s * ldw], 
	    ldw, h[i + cs * ldh], h[i + sn * ldh]);
    *resid = (d__1 = work[i + 1 + s * ldw], fabs(d__1)) / bnrm2;

    if (rank == 0) {
	    printf("iteration no. = %ld, error = %e\n", *iter, *resid);
    }

    if (*resid <= tol) {
	update_(&i, n, &x[1], &h[h_offset], ldh, &work[y * ldw + 1], &
		work[s * ldw + 1], &work[v * ldw + 1], ldw);
	goto L70;
    }
    if (*iter == maxit) {
	goto L50;
    }
    if (i < *restrt) {
	goto L30;
    }

L50:

/*        Compute current solution vector X. */

    update_(restrt, n, &x[1], &h[h_offset], ldh, &work[y * ldw + 1], &
	    work[s * ldw + 1], &work[v * ldw + 1], ldw);

/*        Compute residual vector R, find norm, then check for tolerance. 
*/
/*        (AV is temporary workspace here.) */

    cblas_dcopy(n, &b[1], c__1, &work[av * ldw + 1], c__1);
    (*matvec)(&c_b7, &x[1], &c_b8, &work[av * ldw + 1], particles);
    (*psolve)(&work[r * ldw + 1], &work[av * ldw + 1], particles);
    work[i + 1 + s * ldw] = cblas_dnrm2(n, &work[r * ldw + 1], c__1);
    *resid = work[i + 1 + s * ldw] / bnrm2;
    if (*resid <= tol) {
	goto L70;
    }
    if (*iter == maxit) {
	goto L60;
    }

/*        Restart. */

    goto L10;

L60:

/*     Iteration fails. */

    *info = 1;
    return 0;

L70:

/*     Iteration successful; return. */

    return 0;

/*     End of GMRES */

} /* gmres_ */


/*     =============================================================== */
/* Subroutine */ int update_(long int *i, long int n, double *x, double *h, long int ldh,
                             double *y, double *s, double *v, long int ldv)
{
    /* System generated locals */
    long int h_offset, v_offset;

/*     This routine updates the GMRES iterated solution approximation. */


/*     .. Executable Statements .. */

/*     Solve H*Y = S for upper triangualar H. */

    /* Parameter adjustments */
    v_offset = ldv + 1;
    v -= v_offset;
    --s;
    --y;
    h_offset = ldh + 1;
    h -= h_offset;
    --x;

    /* Function Body */
    cblas_dcopy(*i, &s[1], c__1, &y[1], c__1);
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
                *i, &h[h_offset], ldh, &y[1], c__1);

/*     Compute current solution vector X = X + V*Y. */

    cblas_dgemv(CblasColMajor, CblasNoTrans, n, *i, c_b8, 
                &v[v_offset], ldv, &y[1], c__1, c_b8, &x[1], c__1);

		return 0;

} /* update_ */


/*     ========================================================= */
/* Subroutine */ int basis_(long int *i, long int n, double *h, double *v, long int ldv, double *w)
{
    /* System generated locals */
    long int v_offset, i__1;
    double d__1;

    /* Local variables */
    static long int k;


/*     Construct the I-th column of the upper Hessenberg matrix H */
/*     using the Gram-Schmidt process on V and W. */


    /* Parameter adjustments */
    --w;
    v_offset = ldv + 1;
    v -= v_offset;
    --h;

    /* Function Body */
    i__1 = *i;
    for (k = 1; k <= i__1; ++k) {
	h[k] = cblas_ddot(n, &w[1], c__1, &v[k * ldv + 1], c__1);
	d__1 = -h[k];
	cblas_daxpy(n, d__1, &v[k * ldv + 1], c__1, &w[1], c__1);
/* L10: */
    }
    h[*i + 1] = cblas_dnrm2(n, &w[1], c__1);
    cblas_dcopy(n, &w[1], c__1, &v[(*i + 1) * ldv + 1], c__1);
    d__1 = 1. / h[*i + 1];
    cblas_dscal(n, d__1, &v[(*i + 1) * ldv + 1], c__1);

    return 0;

}
