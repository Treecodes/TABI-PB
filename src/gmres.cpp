#include <iostream>
#include <iomanip>
#include <cmath>

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
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged. The solution is over-written on vector x.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved.
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DROT, DROTG, DSCAL
*  ============================================================
*/

static double dnrm2_(long int n, const double* w);
static void dscal_(long int n, double alpha, double* x);
static double ddot_(long int n, const double* __restrict__ x, const double* __restrict__ y);
static void daxpy_(long int n, double alpha, const double* __restrict__ x, double* __restrict__ y);
static void drot_(double& dx, double& dy, double c, double s);
static void drotg_(double da, double db, double& c, double& s);
static void dtrsv_(long int n, const double* __restrict__ a, long int lda, double* __restrict__ x);
static void dgemv_(long int m, long int n, const double* __restrict__ a, long int lda,
                   const double* __restrict__ x, double* __restrict__ y);

static void update_(long int i, long int n, double* x, const double* h, long int ldh,
                    double* y, const double* s, const double* v, long int ldv);
static void basis_(long int i, long int n, double* h, double* v, long int ldv, double* w);

//*****************************************************************
int Treecode::gmres_(long int n, const double *b, double *x, long int restrt,
                     double* work, long int ldw, double* h, long int ldh,
                     long int& iter, double& resid)
{
    long int maxit = iter;
    double tol = resid;

/*     Store the Givens parameters in matrix H. */
/*     Set initial residual (AV is temporary workspace here). */

    for (long int idx = 0; idx < n; ++idx) work[2 * ldw + idx] = b[idx];

    if (dnrm2_(n, x) != 0.) {
    
        for (long int idx = 0; idx < n; ++idx) work[2 * ldw + idx] = b[idx];
        Treecode::matrix_vector(-1., x, 1., &work[2 * ldw]);
    }

    Treecode::precondition(work, &work[2 * ldw]);

    double bnrm2 = dnrm2_(n, b);
    if (bnrm2 == 0.) bnrm2 = 1.;
    
    if (dnrm2_(n, work) / bnrm2 < tol) {
        return 0;
    }

    iter = 0;

    while (true) {

    /*        Construct the first column of V. */

        for (long int idx = 0; idx < n; ++idx) work[3 * ldw + idx] = work[idx];
        
        double rnorm = dnrm2_(n, &work[3 * ldw]);
        dscal_(n, 1. / rnorm, &work[3 * ldw]);

    /*        Initialize S to the elementary vector E1 scaled by RNORM. */

        work[ldw] = rnorm;
        for (long int k = 1; k < n; ++k) work[k + ldw] = 0.;

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
                drot_(h[k + i * ldh],      h[k + 1 + i * ldh],
                      h[k + restrt * ldh], h[k + (restrt + 1) * ldh]);
            }

        /*           Construct the I-th rotation matrix, and apply it to H so that */
        /*           H(I+1,I) = 0. */
                                  
            drotg_(h[i * (ldh + 1)],      h[i * (ldh + 1) + 1],
                   h[i + (restrt) * ldh], h[i + (restrt + 1) * ldh]);
                             
            drot_ (h[i * (ldh + 1)],      h[i * (ldh + 1) + 1],
                   h[i + (restrt) * ldh], h[i + (restrt + 1) * ldh]);

        /*           Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This */
        /*           gives an approximation of the residual norm. If less than */
        /*           tolerance, update the approximation vector X and quit. */
                            
            drot_(work[i + ldw], work[i + ldw + 1],
                  h[i + (restrt) * ldh], h[i + (restrt + 1) * ldh]);
                            
            resid = std::fabs(work[i + 1 + ldw]) / bnrm2;
            std::cout << "GMRES iteration " << std::setw(3) << iter
                      << ": error = " << std::scientific << resid << std::endl;

            if (resid <= tol) {
                update_(i+1, n, x, h, ldh, &work[2 * ldw], &work[ldw], &work[3 * ldw], ldw);
                return 0;
            }
        }

    /*        Compute current solution vector X. */

        update_(restrt, n, x, h, ldh, &work[2 * ldw], &
                work[ldw], &work[3 * ldw], ldw);

    /*        Compute residual vector R, find norm, then check for tolerance. */

        for (long int idx = 0; idx < n; ++idx) work[2 * ldw + idx] = b[idx];
        
        Treecode::matrix_vector(-1., x, 1., &work[2 * ldw]);
        Treecode::precondition(work, &work[2 * ldw]);
        
        work[restrt + ldw] = dnrm2_(n, work);
        resid = work[restrt + ldw] / bnrm2;
        
        if (resid <= tol) {
            return 0;
        }
        
        if (iter == maxit) {
            return 1;
        }
    } /* Restart. */
}


/*     =============================================================== */
static void update_(long int i, long int n, double* x, const double* h, long int ldh,
                    double* y, const double* s, const double* v, long int ldv)
{
/*     This routine updates the GMRES iterated solution approximation. */
/*     Solve H*Y = S for upper triangualar H. */
/*     Compute current solution vector X = X + V*Y. */

    for (long int idx = 0; idx < i; ++idx) y[idx] = s[idx];
    
    dtrsv_(i, h, ldh, y);
    dgemv_(n, i, v, ldv, y, x);
}


/*     ========================================================= */
static void basis_(long int i, long int n, double* h, double* v, long int ldv, double* w)
{
/*     Construct the I-th column of the upper Hessenberg matrix H */
/*     using the Gram-Schmidt process on V and W. */

    for (long int k = 0; k < i; ++k) {
        h[k] = ddot_(n, w, &v[k * ldv]);
        daxpy_(n, -h[k], &v[k * ldv], w);
    }
    h[i] = dnrm2_(n, w);
    
    for (long int idx = 0; idx < n; ++idx) v[i * ldv + idx] = w[idx];
    dscal_(n, 1. / h[i], &v[i * ldv]);
}


static double dnrm2_(long int n, const double* x)
{
    double norm = 0.;
    for (long int idx = 0; idx < n; ++idx) {
        norm += x[idx] * x[idx];
    }
    return std::sqrt(norm);
}


static void dscal_(long int n, double alpha, double* x)
{
    for (long int idx = 0; idx < n; ++idx) {
        x[idx] *= alpha;
    }
}

static double ddot_(long int n, const double* __restrict__ x,
                    const double* __restrict__ y)
{
    double ddot = 0.;
    for (long int idx = 0; idx < n; ++idx) {
        ddot += x[idx] * y[idx];
    }
    return ddot;
}


static void daxpy_(long int n, double alpha, const double* __restrict__ x,
                   double* __restrict__ y)
{
    for (long int idx = 0; idx < n; ++idx) {
        y[idx] += alpha * x[idx];
    }
}


static void drot_(double& dx, double& dy, double c, double s)
{
/*  applies a plane rotation. */
    double dtemp = c * dx + s * dy;
    dy = c * dy - s * dx;
    dx = dtemp;
}


static void drotg_(double da, double db, double& c, double& s)
{
/*  construct givens plane rotation. */

    double roe = db;
    if (std::abs(da) > std::abs(db)) roe = da;
    double scale = std::abs(da) + std::abs(db);
    
    if (scale != 0.) {
        double d__1 = da / scale;
        double d__2 = db / scale;
        
        double r = scale * std::sqrt(d__1 * d__1 + d__2 * d__2)
                * (roe >= 0. ? 1. : -1.);

        c = da / r;
        s = db / r;
        
    } else {
        c = 1.;
        s = 0.;
    }
}


static void dtrsv_(long int n, const double* __restrict__ a, long int lda,
                   double* __restrict__ x)
{
/*  solve A*x = b, where A is upper triangular */

    for (long int j = n - 1; j >= 0; --j) {
        if (x[j] != 0.) {
            x[j] /= a[j + j*lda];
            double temp = x[j];
            for (long int i = j - 1; i >= 0; --i) {
                x[i] -= temp * a[i + j*lda];
            }
        }
    }
}


static void dgemv_(long int m, long int n, const double* __restrict__ a, long int lda,
                   const double* __restrict__ x, double* __restrict__ y)
{
/*  Form  y = A*x + y */

    for (long int j = 0; j < n; ++j) {
        if (x[j] != 0.) {
            double temp = x[j];
            for (long int i = 0; i < m; ++i) {
                y[i] += temp * a[i + j*lda];
            }
        }
    }
}
