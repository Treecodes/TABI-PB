/**************************************************************************
* FILE NAME: treecode.c                                                   *
*                                                                         *
* PURPOSE: Contains all treecode-related functions and variables,         *
*          including treecode initialization and finalization functions   *
*          that interface with tabipb.c, and matrix-vector multiplication *
*          and solve functions that interface with run_gmres.c            *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//#include <Accelerate/Accelerate.h>
//#include <lapacke.h>
#include <mkl_lapacke.h>

#ifdef MPI_ENABLED
    #include <mpi.h>
#endif

#include "treecode_tabipb_interface.h"
#include "treecode_gmres_interface.h"
#include "utilities.h"
#include "partition.h"

#include "global_params.h"
#include "array.h"
#include "tree_node_struct.h"
#include "particle_struct.h"
#include "TABIPBstruct.h"


/* runtime treecode parameters */
static int s_numpars;
static int s_order;
static int s_max_per_leaf;
static double theta;

/* runtime physical parameters */
static double s_kappa;
static double s_kappa2;
static double s_eps;

/* variables for tracking tree information */
static int s_min_level;
static int s_max_level;

/* global variables for taylor expansions */
static int s_torder_flat;

static double *s_cf1 = NULL;
static double *s_cf2 = NULL;
static double *s_cf3 = NULL;

static double ***s_a = NULL;
static double ***s_b = NULL;

static double ****s_der_coeff = NULL;
static int s_kk[16][3];

/* these point to arrays located in TreeParticles */
static double **s_particle_position = NULL;
static double **s_particle_normal = NULL;
static double *s_particle_area = NULL;
static double *s_source_term = NULL;

/* global variables used when computing potential/force */
static double s_target_position[3];
static double s_target_normal[3];

static double ***s_target_charge = NULL;
static double ***s_source_charge = NULL;

/* global variables for reordering arrays */
static int *s_order_arr = NULL;

/* root node of tree */
static TreeNode *s_tree_root = NULL;


/* internal functions */
static int s_Setup(double xyz_limits[6]);
static int s_CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6],
                        int level);
static int s_PartitionEight(double xyzmms[6][8], double xl, double yl,
                            double zl, double lmax, double x_mid, double y_mid,
                            double z_mid, int ind[8][2]);
static int s_ComputePBKernel(double *phi);
static int s_ComputeAllMoments(TreeNode *p, int ifirst);
static int s_ComputeMoments(TreeNode *p);
static int s_RunTreecode(TreeNode *p, double *tpoten_old,
                         double tempq[2][16], double peng[2]);
static int s_ComputeTreePB(TreeNode *p, double tempq[2][16], double peng[2]);
static int s_ComputeCoeffs(TreeNode *p);
static int s_ComputeCoeffsCoulomb(TreeNode *p);
static int s_ComputeDirectPB(int ibeg, int iend, double *tpoten_old,
                             double peng[2]);
static int s_RemoveMoments(TreeNode *p);
static int s_RemoveNode(TreeNode *p);

/* internal preconditioning functions */
static void leaflength(TreeNode *p, int idx, int *nrow);
static int lu_decomp(double **A, int N, int *ipiv);
static void lu_solve(double **matrixA, int N, int *ipiv, double *rhs);



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* TreecodeInitialization and Finalization are used by       * * * */
/* tabipb() to interface with the treecode                   * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
int TreecodeInitialization(TABIPBparm *parm, TreeParticles *particles)
{
    /* set up variables used in treecode */
    /* local variables*/
    int level, i, j, k, mm, nn, idx, ijk[3], ierr;

    /* variables needed for reorder */
    double *temp_area, *temp_source;
    double **temp_normal;
    
    double xyz_limits[6];
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

    if (rank == 0) {
        printf("\nInitializing treecode...\n");
    }

    /* setting variables global to file */
    s_numpars = particles->num_particles;
    s_order = parm->order;
    s_max_per_leaf = parm->maxparnode;
    theta = parm->theta;
    
    s_kappa = parm->kappa;
    s_kappa2 = parm->kappa2;
    s_eps = parm->eps;
    
    s_torder_flat = (s_order+1) * (s_order+2) * (s_order+3) / 6;
    
    s_min_level = 50000;
    s_max_level = 0;

    level = 0;

    /* allocate s_der_coeff */
    
    make_4array(s_der_coeff, s_order+1, s_order+1, s_order+1, 16);

    /* initialize s_kk */
    s_kk[0][0]=0; s_kk[0][1]=0; s_kk[0][2]=0; /* Original Kernel */

    s_kk[1][0]=1; s_kk[1][1]=0; s_kk[1][2]=0; /* 1st Order */
    s_kk[2][0]=0; s_kk[2][1]=1; s_kk[2][2]=0;
    s_kk[3][0]=0; s_kk[3][1]=0; s_kk[3][2]=1;

    s_kk[4][0]=1; s_kk[4][1]=0; s_kk[4][2]=0;
    s_kk[5][0]=0; s_kk[5][1]=1; s_kk[5][2]=0;
    s_kk[6][0]=0; s_kk[6][1]=0; s_kk[6][2]=1;

    s_kk[7][0]=2; s_kk[7][1]=0; s_kk[7][2]=0;
    s_kk[8][0]=1; s_kk[8][1]=1; s_kk[8][2]=0;
    s_kk[9][0]=1; s_kk[9][1]=0; s_kk[9][2]=1;
    s_kk[10][0]=1; s_kk[10][1]=1; s_kk[10][2]=0;
    s_kk[11][0]=0; s_kk[11][1]=2; s_kk[11][2]=0;
    s_kk[12][0]=0; s_kk[12][1]=1; s_kk[12][2]=1;
    s_kk[13][0]=1; s_kk[13][1]=0; s_kk[13][2]=1;
    s_kk[14][0]=0; s_kk[14][1]=1; s_kk[14][2]=1;
    s_kk[15][0]=0; s_kk[15][1]=0; s_kk[15][2]=2;

    /* the adjustment of s_der_coeff for the recurrence relation */

    for (i = 0; i < s_order+1; i++) {
        for (j = 0; j < s_order+1; j++) {
            for (k = 0; k < s_order+1; k++) {
                for (idx = 0; idx < 16; idx++) {
                    s_der_coeff[i][j][k][idx] = 1.0;
                }
            }
        }
    }

    for (i = 0; i < s_order+1; i++) {
        ijk[0] = i;
        for (j = 0; j < s_order+1; j++) {
            ijk[1] = j;
            for (k = 0; k < s_order+1; k++) {
                ijk[2] = k;
                for (idx = 0; idx < 16; idx++) {
                    for (mm = 0; mm < 3; mm++) {
                        if (s_kk[idx][mm] != 0) {
                            for (nn = 0; nn < s_kk[idx][mm]; nn++)
                                s_der_coeff[i][j][k][idx] *= (ijk[mm] + (nn+1));
                        }
                    }
                }
            }   
        }
    }


    for (i = 0; i < s_order+1; i++) {
        for (j = 0; j < s_order+1; j++) {
            for (k = 0; k < s_order+1; k++) {
                for (idx = 0; idx < 16; idx++) {
                    s_der_coeff[i][j][k][idx] *= ONE_OVER_4PI;
                }
            }
        }
    }
    
    s_particle_position = particles->position;
    s_particle_normal = particles->normal;
    s_particle_area = particles->area;
    s_source_term = particles->source_term;

    make_matrix(temp_normal, 3, s_numpars);
    make_vector(temp_area, s_numpars);
    make_vector(temp_source, 2 * s_numpars);
    

/* Call SETUP to allocate arrays for Taylor expansions */
/* and setup global variables. Also, copy variables into global copy arrays. */
    s_Setup(xyz_limits);

    s_tree_root = (TreeNode*)calloc(1, sizeof(TreeNode));

    s_CreateTree(s_tree_root, 0, s_numpars-1, xyz_limits, level);
    
    if (rank == 0) {
        printf("Created tree for %d particles with max %d per node.\n\n",
               s_numpars, s_max_per_leaf);
    }

    memcpy(temp_normal[0], s_particle_normal[0], s_numpars*sizeof(double));
    memcpy(temp_normal[1], s_particle_normal[1], s_numpars*sizeof(double));
    memcpy(temp_normal[2], s_particle_normal[2], s_numpars*sizeof(double));
    memcpy(temp_area, s_particle_area, s_numpars*sizeof(double));
    memcpy(temp_source, s_source_term, 2*s_numpars*sizeof(double));
    
    for (i = 0; i < s_numpars; i++) {
        s_particle_normal[0][i]    = temp_normal[0][s_order_arr[i]];
        s_particle_normal[1][i]    = temp_normal[1][s_order_arr[i]];
        s_particle_normal[2][i]    = temp_normal[2][s_order_arr[i]];
        s_particle_area[i]         = temp_area[s_order_arr[i]];
        s_source_term[i]           = temp_source[s_order_arr[i]];
        s_source_term[i + s_numpars] = temp_source[s_order_arr[i] + s_numpars];
    }

    free_matrix(temp_normal);
    free_vector(temp_area);
    free_vector(temp_source);

    make_3array(s_target_charge, s_numpars, 2, 16);
    make_3array(s_source_charge, s_numpars, 2, 16);

    return 0;
}
/**********************************************************/


/********************************************************/
int TreecodeFinalization(TreeParticles *particles)
{

    int i, ierr;
    double *temp_area, *temp_source, *temp_xvct;
    double **temp_normal, **temp_position;
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

/***********reorder particles*************/

    make_matrix(temp_position, 3, s_numpars);
    make_matrix(temp_normal, 3, s_numpars);
    make_vector(temp_area, s_numpars);
    make_vector(temp_source, 2 * s_numpars);
    make_vector(temp_xvct, 2 * s_numpars);

    memcpy(temp_position[0], particles->position[0], s_numpars*sizeof(double));
    memcpy(temp_position[1], particles->position[1], s_numpars*sizeof(double));
    memcpy(temp_position[2], particles->position[2], s_numpars*sizeof(double));
    memcpy(temp_normal[0], particles->normal[0], s_numpars*sizeof(double));
    memcpy(temp_normal[1], particles->normal[1], s_numpars*sizeof(double));
    memcpy(temp_normal[2], particles->normal[2], s_numpars*sizeof(double));
    memcpy(temp_area, particles->area, s_numpars*sizeof(double));
    memcpy(temp_source, particles->source_term, 2*s_numpars*sizeof(double));
    memcpy(temp_xvct, particles->xvct, 2*s_numpars*sizeof(double));
    
    for (i = 0; i < s_numpars; i++) {
        particles->position[0][s_order_arr[i]]  = temp_position[0][i];
        particles->position[1][s_order_arr[i]]  = temp_position[1][i];
        particles->position[2][s_order_arr[i]]  = temp_position[2][i];
        particles->normal[0][s_order_arr[i]]    = temp_normal[0][i];
        particles->normal[1][s_order_arr[i]]    = temp_normal[1][i];
        particles->normal[2][s_order_arr[i]]    = temp_normal[2][i];
        particles->area[s_order_arr[i]]         = temp_area[i];
        particles->source_term[s_order_arr[i]]  = temp_source[i];
        particles->source_term[s_order_arr[i] + s_numpars]
                                                = temp_source[i + s_numpars];
        particles->xvct[s_order_arr[i]]         = temp_xvct[i];
        particles->xvct[s_order_arr[i] + s_numpars]
                                                = temp_xvct[i + s_numpars];
    }

    free_matrix(temp_position);
    free_matrix(temp_normal);
    free_vector(temp_area);
    free_vector(temp_source);
    free_vector(temp_xvct);

/***********treecode_initialization*******/
    free_4array(s_der_coeff);

    free_3array(s_target_charge);
    free_3array(s_source_charge);
    
/***********clean tree structure**********/

    s_RemoveNode(s_tree_root);
    free(s_tree_root);

/***********variables in setup************/
    free_vector(s_cf1);
    free_vector(s_cf2);
    free_vector(s_cf3);
    
    free_vector(s_a);
    free_vector(s_b);

    free_vector(s_order_arr);
/*****************************************/

    if (rank == 0) {
        printf("\nTABIPB tree structure has been deallocated.\n\n");
    }

    return 0;
}
/**********************************************************/



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* matvec and psolve are the functions called by GMRes             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
int matvec(double *alpha, double *tpoten_old, double *beta, double *tpoten)
{
/* the main part of treecode */
/* in gmres *matvec(Alpha, X, Beta, Y) where y := alpha*A*x + beta*y */
  /* local variables */
    int i, j, k, ii, ierr;
    double temp_x, temp_area;
    double temp_charge[2][16];
    double pre1, pre2;
    double peng[2], peng_old[2];
    double *tpoten_temp;
    
    int particles_per_process;
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif
    
    make_vector(tpoten_temp, 2 * s_numpars);
    memcpy(tpoten_temp, tpoten, 2 * s_numpars * sizeof(double));
    memset(tpoten, 0, 2 * s_numpars * sizeof(double));
    
    s_ComputePBKernel(tpoten_old);
    
  /* Generate the moments if not allocated yet */
    s_ComputeAllMoments(s_tree_root, 1);

    pre1 = 0.50 * (1.0 + s_eps);
    pre2 = 0.50 * (1.0 + 1.0/s_eps);
    
    particles_per_process = s_numpars / num_procs;

    for (ii = 0; ii <= particles_per_process; ii++) {
        i = ii * num_procs + rank;
        
        if (i < s_numpars) {
            peng[0] = 0.0;
            peng[1] = 0.0;
            peng_old[0] = tpoten_old[i];
            peng_old[1] = tpoten_old[i+s_numpars];
            s_target_position[0] = s_particle_position[0][i];
            s_target_position[1] = s_particle_position[1][i];
            s_target_position[2] = s_particle_position[2][i];
            s_target_normal[0] = s_particle_normal[0][i];
            s_target_normal[1] = s_particle_normal[1][i];
            s_target_normal[2] = s_particle_normal[2][i];
        
            for (j = 0; j < 2; j++) {
                for (k = 0; k < 16; k++) {
                    temp_charge[j][k] = s_target_charge[i][j][k];
                }
            }

      /* remove the singularity */
            temp_x = s_particle_position[0][i];
            temp_area = s_particle_area[i];
            s_particle_position[0][i] += 100.123456789;
            s_particle_area[i] = 0.0;

      /* start to use Treecode */
            s_RunTreecode(s_tree_root, tpoten_old, temp_charge, peng);

            tpoten[i] = tpoten_temp[i] * *beta
                      + (pre1 * peng_old[0] - peng[0]) * *alpha;
            tpoten[s_numpars+i] = tpoten_temp[s_numpars+i] * *beta
                                + (pre2 * peng_old[1] - peng[1]) * *alpha;

            s_particle_position[0][i] = temp_x;
            s_particle_area[i] = temp_area;
        }
    }
    
#ifdef MPI_ENABLED
    ierr = MPI_Allreduce(MPI_IN_PLACE, tpoten, 2 * s_numpars,
                         MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    free_vector(tpoten_temp);

    s_RemoveMoments(s_tree_root);

    return 0;
}
/**********************************************************/


/**********************************************************/
int psolve(double *z, double *r)
{
/* r as original while z as scaled */

        int i;
        double scale1, scale2;
        scale1 = 0.5 * (1.0 + s_eps);
        scale2 = 0.5 * (1.0 + 1.0/s_eps);

        for (i = 0; i < s_numpars; i++) {
                z[i] = r[i]/scale1;
                z[i + s_numpars] = r[i + s_numpars]/scale2;
        }

        return 0;
}
/**********************************************************/


/**********************************************************/
int psolve_precond(double *z, double *r)
{
/* r as original while z as scaled */

    int i, j, idx = 0, nrow = 0, nrow2, ibeg = 0, iend = 0;
    int *ipiv, inc, nrhs, info;
    //double **matrixA;
    double *columnMajorA, *rhs;
    double L1, L2, L3, L4, area;
    double tp[3], tq[3], sp[3], sq[3];
    double r_s[3], rs, irs, sumrs;
    double G0, kappa_rs, exp_kappa_rs, Gk;
    double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
    double G10, G20, G1, G2, G3, G4;
    double pre1, pre2;

    pre1 = 0.5*(1.0+s_eps);
    pre2 = 0.5*(1.0+1.0/s_eps);
  
    //make_matrix(matrixA, 2*s_max_per_leaf, 2*s_max_per_leaf);
    make_vector(columnMajorA, 4*s_max_per_leaf*s_max_per_leaf);
    make_vector(ipiv, 2*s_max_per_leaf);
    make_vector(rhs, 2*s_max_per_leaf);

    while (idx < s_numpars) {
        leaflength(s_tree_root, idx, &nrow);
        nrow2 = nrow*2;
        ibeg  = idx;
        iend  = idx + nrow - 1;

        memset(columnMajorA, 0, nrow2*nrow2*sizeof(double));
        memset(ipiv, 0, nrow2*sizeof(int));
        memset(rhs, 0, nrow2*sizeof(double));

        for (i = ibeg; i <= iend; i++) {
            tp[0] = s_particle_position[0][i];
            tp[1] = s_particle_position[1][i];
            tp[2] = s_particle_position[2][i];
            tq[0] = s_particle_normal[0][i];
            tq[1] = s_particle_normal[1][i];
            tq[2] = s_particle_normal[2][i];

            for (j = ibeg; j < i; j++) {
                sp[0] = s_particle_position[0][j];
                sp[1] = s_particle_position[1][j];
                sp[2] = s_particle_position[2][j];
                sq[0] = s_particle_normal[0][j];
                sq[1] = s_particle_normal[1][j];
                sq[2] = s_particle_normal[2][j];

                r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
                sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
                rs = sqrt(sumrs);
                irs = 1.0/rs;
                G0 = ONE_OVER_4PI * irs;
                kappa_rs = s_kappa * rs;
                exp_kappa_rs = exp(-kappa_rs);
                Gk = exp_kappa_rs * G0;

                cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
                cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
                tp1 = G0* irs;
                tp2 = (1.0 + kappa_rs) * exp_kappa_rs;

                G10 = cos_theta0 * tp1;
                G20 = tp2 * G10;

                G1 = cos_theta * tp1;
                G2 = tp2 * G1;

                dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
                G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
                G4 = tp2*G3 - s_kappa2*cos_theta0*cos_theta*Gk;

                area = s_particle_area[j];

                L1 = G1 - s_eps*G2;
                L2 = G0 - Gk;
                L3 = G4 - G3;
                L4 = G10 - G20/s_eps;

                //matrixA[i-ibeg][j-ibeg] = -L1*area;
                //matrixA[i-ibeg][j+nrow-ibeg] = -L2*area;
                //matrixA[i+nrow-ibeg][j-ibeg] = -L3*area;
                //matrixA[i+nrow-ibeg][j+nrow-ibeg] = -L4*area;
                columnMajorA[(j-ibeg)*nrow2 + i-ibeg] = -L1*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i-ibeg] = -L2*area;
                columnMajorA[(j-ibeg)*nrow2 + i+nrow-ibeg] = -L3*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i+nrow-ibeg] = -L4*area;
            }

            //matrixA[i-ibeg][i-ibeg] = pre1;
            //matrixA[i+nrow-ibeg][i+nrow-ibeg] = pre2;
            columnMajorA[(i-ibeg)*nrow2 + i-ibeg] = pre1;
            columnMajorA[(i+nrow-ibeg)*nrow2 + i+nrow-ibeg] = pre2;

            for (j = i+1; j <= iend; j++) {
                sp[0] = s_particle_position[0][j];
                sp[1] = s_particle_position[1][j];
                sp[2] = s_particle_position[2][j];
                sq[0] = s_particle_normal[0][j];
                sq[1] = s_particle_normal[1][j];
                sq[2] = s_particle_normal[2][j];

                r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
                sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
                rs = sqrt(sumrs);
                irs = 1.0/rs;
                G0 = ONE_OVER_4PI * irs;
                kappa_rs = s_kappa * rs;
                exp_kappa_rs = exp(-kappa_rs);
                Gk = exp_kappa_rs * G0;

                cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
                cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
                tp1 = G0 * irs;
                tp2 = (1.0 + kappa_rs) * exp_kappa_rs;

                G10 = cos_theta0 * tp1;
                G20 = tp2 * G10;

                G1 = cos_theta * tp1;
                G2 = tp2 * G1;

                dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
                G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
                G4 = tp2*G3 - s_kappa2*cos_theta0*cos_theta*Gk;

                area = s_particle_area[j];

                L1 = G1 - s_eps*G2;
                L2 = G0 - Gk;
                L3 = G4 - G3;
                L4 = G10 - G20/s_eps;

                //matrixA[i-ibeg][j-ibeg] = -L1*area;
                //matrixA[i-ibeg][j+nrow-ibeg] = -L2*area;
                //matrixA[i+nrow-ibeg][j-ibeg] = -L3*area;
                //matrixA[i+nrow-ibeg][j+nrow-ibeg] = -L4*area;
                columnMajorA[(j-ibeg)*nrow2 + i-ibeg] = -L1*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i-ibeg] = -L2*area;
                columnMajorA[(j-ibeg)*nrow2 + i+nrow-ibeg] = -L3*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i+nrow-ibeg] = -L4*area;
            }
        }

        for (i = 0; i < nrow; i++) {
            rhs[i] = r[i+ibeg];
            rhs[i+nrow] = r[i+ibeg+s_numpars];
        }

        // Jiahui's implementation of LU decomposition
        //inc = lu_decomp(matrixA, nrow2, ipiv);
        //lu_solve(matrixA, nrow2, ipiv, rhs);

        // Apple Accelerate implementation of LAPACK LU decomposition
        //nrhs = 1;
        //dgesv_(&nrow2, &nrhs, columnMajorA, &nrow2, ipiv, rhs, &nrow2, &info);

        // LAPACKE implementation of LAPACK LU decomposition
        LAPACKE_dgesv(LAPACK_COL_MAJOR, nrow2, 1, columnMajorA, nrow2, ipiv, rhs, nrow2);

        for (i = 0; i < nrow; i++) {
            z[i+ibeg] = rhs[i];
            z[i+ibeg+s_numpars] = rhs[i+nrow];
        }

        idx += nrow;
    }

    //free_matrix(matrixA);
    free_vector(columnMajorA);
    free_vector(rhs);
    free_vector(ipiv);
  
    return 0;
}
/**********************************************************/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* PRIVATE internal treecode functions                       * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
static int s_Setup(double xyz_limits[6])
{
/* SETUP allocates and initializes arrays needed for the Taylor expansion.
 Also, global variables are set and the Cartesian coordinates of
 the smallest box containing the particles is determined. The particle
 positions and charges are copied so that they can be restored upon exit.*/
 
    int i, j, k, ierr;
    double t1;

    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif
    
    if (rank == 0) {
        printf("Setting up arrays for Taylor expansion...\n");
    }
    
/* allocate global Taylor expansion variables */
  
    make_vector(s_cf1, s_order+2);
    make_vector(s_cf2, s_order+2);
    make_vector(s_cf3, s_order+2);
    
    make_3array(s_a, s_order+5, s_order+5, s_order+5);
    make_3array(s_b, s_order+5, s_order+5, s_order+5);

    for (i = 0; i < s_order+5; i++) {
        for (j = 0; j < s_order+5; j++) {
            for (k = 0; k < s_order+5; k++) {
                s_a[i][j][k] = 0.0;
                s_b[i][j][k] = 0.0;
            }
        }
    }

    for (i = 0; i < s_order+2; i++) {
        t1 = 1.0 / (i+1.0);
        s_cf1[i] = t1;
        s_cf2[i] = 1.0 - 0.5*t1;
        s_cf3[i] = 1.0 - t1;
    }

/* find bounds of Cartesion box enclosing the particles */

    xyz_limits[0] = MinVal(s_particle_position[0], s_numpars);
    xyz_limits[1] = MaxVal(s_particle_position[0], s_numpars);
    xyz_limits[2] = MinVal(s_particle_position[1], s_numpars);
    xyz_limits[3] = MaxVal(s_particle_position[1], s_numpars);
    xyz_limits[4] = MinVal(s_particle_position[2], s_numpars);
    xyz_limits[5] = MaxVal(s_particle_position[2], s_numpars);

    //printf("x-limits of box: %f, %f\n", xyz_limits[0], xyz_limits[1]);
    //printf("y-limits of box: %f, %f\n", xyz_limits[2], xyz_limits[3]);
    //printf("z-limits of box: %f, %f\n", xyz_limits[4], xyz_limits[5]);

    make_vector(s_order_arr, s_numpars);

    for (i = 0; i < s_numpars; i++) {
        s_order_arr[i] = i;
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6],
                        int level)
{
/*CREATE_TREE recursively create the tree structure. Node P is
  input, which contains particles indexed from IBEG to IEND. After
  the node parameters are set subdivision occurs if IEND-IBEG+1 > s_max_per_leaf.
  Real array XYZMM contains the min and max values of the coordinates
  of the particle in P, thus defining the box. */

  /* local variables */
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int ind[8][2];
    double xyzmms[6][8];
    int i, j, loclev, numposchild;
    double lxyzmm[6];

/* set node fields: number of particles, exist_ms and xyz bounds */
    p->numpar = iend-ibeg+1;
    p->exist_ms = 0;

    p->x_min = xyzmm[0];
    p->x_max = xyzmm[1];
    p->y_min = xyzmm[2];
    p->y_max = xyzmm[3];
    p->z_min = xyzmm[4];
    p->z_max = xyzmm[5];
    
/* compute aspect ratio */
    xl = p->x_max-p->x_min;
    yl = p->y_max-p->y_min;
    zl = p->z_max-p->z_min;

    lmax = xl;
    if (lmax < yl) lmax = yl;
    if (lmax < zl) lmax = zl;

    t1 = lmax;
    t2 = xl;
    if (t2 > yl) t2 = yl;
    if (t2 > zl) t2 = zl;

    if (t2 != 0.0) {
        p->aspect = t1/t2;
    } else {
        p->aspect = 0.0;
    }
    
/* midpoint coordinates, RADIUS and SQRADIUS */
    p->x_mid = (p->x_max + p->x_min) / 2.0;
    p->y_mid = (p->y_max + p->y_min) / 2.0;
    p->z_mid = (p->z_max + p->z_min) / 2.0;
    t1 = p->x_max - p->x_mid;
    t2 = p->y_max - p->y_mid;
    t3 = p->z_max - p->z_mid;
    p->radius = sqrt(t1*t1 + t2*t2 + t3*t3);

/* set particle limits, tree level of node, and nullify children pointers */
    p->ibeg = ibeg;
    p->iend = iend;
    p->level = level;
    if (s_max_level < level) s_max_level = level;

    p->num_children = 0;

    make_vector(p->child, 8);
    for (i = 0; i < 8; i++) {
        p->child[i] = (TreeNode*)calloc(1, sizeof(TreeNode));
    }
    

    if (p->numpar > s_max_per_leaf) {
/* set IND array to 0 and then call PARTITION routine. IND array holds indices
 * of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1 */

        xyzmms[0][0] = p->x_min;
        xyzmms[1][0] = p->x_max;
        xyzmms[2][0] = p->y_min;
        xyzmms[3][0] = p->y_max;
        xyzmms[4][0] = p->z_min;
        xyzmms[5][0] = p->z_max;
        
        for (i = 0; i < 8; i++) {
            ind[i][0] = 0;
            ind[i][1] = 0;
        }
        
        ind[0][0] = ibeg;
        ind[0][1] = iend;
        x_mid = p->x_mid;
        y_mid = p->y_mid;
        z_mid = p->z_mid;

        numposchild = s_PartitionEight(xyzmms, xl, yl, zl, lmax,
                                       x_mid, y_mid, z_mid, ind);
        
/* Shrink the box */
        for (i = 0; i < 8; i++) {
            if (ind[i][0] < ind[i][1]) {
                xyzmms[0][i] = MinVal(&s_particle_position[0][ind[i][0]],
                                      ind[i][1]-ind[i][0]);
                xyzmms[1][i] = MaxVal(&s_particle_position[0][ind[i][0]],
                                      ind[i][1]-ind[i][0]);
                xyzmms[2][i] = MinVal(&s_particle_position[1][ind[i][0]],
                                      ind[i][1]-ind[i][0]);
                xyzmms[3][i] = MaxVal(&s_particle_position[1][ind[i][0]],
                                      ind[i][1]-ind[i][0]);
                xyzmms[4][i] = MinVal(&s_particle_position[2][ind[i][0]],
                                      ind[i][1]-ind[i][0]);
                xyzmms[5][i] = MaxVal(&s_particle_position[2][ind[i][0]],
                                      ind[i][1]-ind[i][0]);
            }
        }
/* create children if indicated and store info in parent */
        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                p->num_children = p->num_children + 1;
                for (j = 0; j < 6; j++) {
                    lxyzmm[j] = xyzmms[j][i];
                }
                s_CreateTree(p->child[p->num_children-1],
                             ind[i][0], ind[i][1], lxyzmm, loclev);
            }
        }
    } else {
        if (level < s_min_level) {
            s_min_level = level;
        }
    }

    return 0;
}
/**********************************************************/


/********************************************************/
static int s_PartitionEight(double xyzmms[6][8], double xl, double yl,
                            double zl, double lmax, double x_mid, double y_mid,
                            double z_mid, int ind[8][2])
{
/* PARTITION_8 determines the particle indices of the eight sub boxes
 * containing the particles after the box defined by particles I_BEG
 * to I_END is divided by its midpoints in each coordinate direction.
 * The determination of the indices is accomplished by the subroutine
 * PARTITION. A box is divided in a coordinate direction as long as the
 * resulting aspect ratio is not too large. This avoids the creation of
 * "narrow" boxes in which Talyor expansions may become inefficient.
 * On exit the INTEGER array IND (dimension 8 x 2) contains
 * the indice limits of each new box (node) and NUMPOSCHILD the number
 * of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
 * that box J is empty.*/
    int temp_ind, i, j;
    double critlen;
    int numposchild;

    numposchild = 1;
    critlen = lmax/sqrt(2.0);

    if (xl >= critlen) {
        temp_ind = Partition(s_particle_position[0], s_particle_position[1],
                             s_particle_position[2], s_order_arr,
                             ind[0][0], ind[0][1], x_mid);
        ind[1][0] = temp_ind+1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;
        for (i = 0; i < 6; i++) {
            xyzmms[i][1] = xyzmms[i][0];
        }
        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        numposchild *= 2;
    }

    if (yl >= critlen) {
        for (i = 0; i < numposchild; i++) {
            temp_ind = Partition(s_particle_position[1], s_particle_position[0],
                                 s_particle_position[2], s_order_arr,
                                 ind[i][0], ind[i][1], y_mid);
            ind[numposchild+i][0] = temp_ind+1;
            ind[numposchild+i][1] = ind[i][1];
            ind[i][1] = temp_ind;
            for (j = 0; j < 6; j++) {
                xyzmms[j][numposchild+i] = xyzmms[j][i];
            }
            xyzmms[3][i] = y_mid;
            xyzmms[2][numposchild+i] = y_mid;
        }
        numposchild *= 2;
    }

    if (zl >= critlen) {
        for (i = 0; i < numposchild; i++) {
            temp_ind = Partition(s_particle_position[2], s_particle_position[0],
                                 s_particle_position[1], s_order_arr,
                                 ind[i][0], ind[i][1], z_mid);
            ind[numposchild+i][0] = temp_ind+1;
            ind[numposchild+i][1] = ind[i][1];
            ind[i][1] = temp_ind;
            for (j = 0; j < 6; j++) {
                xyzmms[j][numposchild+i] = xyzmms[j][i];
            }
            xyzmms[5][i] = z_mid;
            xyzmms[4][numposchild+i] = z_mid;
        }
        numposchild *= 2;
    }

    return (numposchild);
}
/********************************************************/


/********************************************************/
static int s_ComputePBKernel(double *phi)
{

    int i, ikp, iknl, ixyz, jxyz, indx;

    for (i = 0; i < s_numpars; i++) {
        indx = 0;
        for (ikp = 0; ikp < 2; ikp++) {
            s_target_charge[i][ikp][indx] = 1.0;
            s_source_charge[i][ikp][indx] = s_particle_area[i]
                                          * phi[s_numpars+i];
        }

        for (iknl = 0; iknl < 2; iknl++) {
            for (ixyz = 0; ixyz < 3; ixyz++) {
                indx += 1;
                for (ikp = 0; ikp < 2; ikp++) {
                    s_target_charge[i][ikp][indx] = 1.0 * (1-iknl)
                                                  + s_particle_normal[ixyz][i]
                                                  * iknl;
                    s_source_charge[i][ikp][indx] = (s_particle_normal[ixyz][i]
                                                  * (1-iknl) + 1.0 * iknl)
                                                  * s_particle_area[i]
                                                  * phi[iknl*s_numpars+i];
                }
            }
        }

        for (ixyz = 0; ixyz < 3; ixyz++) {
            for (jxyz = 0; jxyz < 3; jxyz++) {
                indx += 1;
                for (ikp = 0; ikp < 2; ikp++) {
                    s_target_charge[i][ikp][indx] =  s_particle_normal[jxyz][i];
                    s_source_charge[i][ikp][indx] = -s_particle_normal[ixyz][i]
                                                  * s_particle_area[i] * phi[i];
                }
            }
        }
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeAllMoments(TreeNode *p, int ifirst)
{
/* REMOVE_NODE recursively removes each node from the tree and deallocates
 * its memory for MS array if it exits. */

    int i;
    
    if (p->exist_ms == 0 && ifirst == 0) {
        make_matrix(p->ms, 16, s_torder_flat);
        s_ComputeMoments(p);
        p->exist_ms = 1;
    }

    if (p->num_children > 0) {
        for (i = 0; i < p->num_children; i++) {
            s_ComputeAllMoments(p->child[i], 0);
        }
    }

  return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeMoments(TreeNode *p)
{
/* COMP_MS computes the moments for node P needed in the Taylor
 * approximation */

    int ii, i, j, k1, k2, k3, n;
    double dx, dy, dz, tx, ty, tz, txyz;
    
    for (n = 0; n < 16; n++) {
        for (i = 0; i < s_torder_flat; i++) {
            p->ms[n][i] = 0.0;
        }
    }

    for (i = p->ibeg; i < p->iend+1; i++) {
        dx = s_particle_position[0][i]-p->x_mid;
        dy = s_particle_position[1][i]-p->y_mid;
        dz = s_particle_position[2][i]-p->z_mid;
        for (j = 0; j < 7; j++) {
        
            ii = 0;
            
            tx = 1.0;
            for (k1 = 0; k1 < s_order+1; k1++) {
                ty=1.0;
                for (k2 = 0; k2 < s_order+1-k1; k2++) {
                    tz=1.0;
                    for (k3 = 0; k3 < s_order+1-k1-k2; k3++) {

                        txyz = tx*ty*tz;
                        p->ms[j][ii] += s_source_charge[i][0][j] * txyz;
                        ii++;
                        
                        tz = tz*dz;
                    }
                    ty = ty*dy;
                }
                tx = tx*dx;
            }
        }

        ii = 0;

        tx = 1.0;
        for (k1 = 0; k1 < s_order+1; k1++) {
            ty = 1.0;
            for (k2 = 0; k2 < s_order+1-k1; k2++) {
                tz = 1.0;
                for (k3 = 0; k3 < s_order+1-k1-k2; k3++) {

                    txyz = tx*ty*tz;
                    p->ms[7][ii] += s_source_charge[i][0][7] * txyz;
                    p->ms[10][ii] += s_source_charge[i][0][10] * txyz;
                    p->ms[13][ii] += s_source_charge[i][0][13] * txyz;
                    ii++;

                    tz *= dz;
                }
                ty *= dy;
            }
            tx *= dx;
        }
    }
    
    ii = 0;

    for (k1 = 0; k1 < s_order+1; k1++) {
        for (k2 = 0; k2 < s_order+1-k1; k2++) {
            for (k3 = 0; k3 < s_order+1-k1-k2; k3++) {
                p->ms[8][ii] = p->ms[7][ii];
                p->ms[9][ii] = p->ms[7][ii];
                p->ms[11][ii] = p->ms[10][ii];
                p->ms[12][ii] = p->ms[10][ii];
                p->ms[14][ii] = p->ms[13][ii];
                p->ms[15][ii] = p->ms[13][ii];
                ii++;
            }
        }
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_RunTreecode(TreeNode *p, double *tpoten_old, double tempq[2][16],
                         double peng[2])
{
  /* RunTreecode() is self recurrence function */
  
    double tx, ty, tz, dist, pengchild[2];
    int i;


  /* determine DISTSQ for MAC test */
    tx = p->x_mid - s_target_position[0];
    ty = p->y_mid - s_target_position[1];
    tz = p->z_mid - s_target_position[2];
    dist = sqrt(tx*tx + ty*ty + tz*tz);

  /* initialize potential energy */
    peng[0] = 0.0; 
    peng[1] = 0.0;

/* If MAC is accepted and there is more than 1 particale in the */
/* box use the expansion for the approximation. */

    if (p->radius < dist*theta && p->numpar > 40) {
        s_ComputeTreePB(p, tempq, peng);
    } else {
        if (p->num_children == 0) {
            s_ComputeDirectPB(p->ibeg, p->iend, tpoten_old, peng);
        } else {
      /* If MAC fails check to see if there are children. If not, perform */
      /* direct calculation.  If there are children, call routine */
      /* recursively for each. */
            for (i = 0; i < p->num_children; i++) {
                pengchild[0] = 0.0; 
                pengchild[1] = 0.0;
                s_RunTreecode(p->child[i], tpoten_old, tempq, pengchild);
                peng[0] += pengchild[0];
                peng[1] += pengchild[1];
            }
        }
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeTreePB(TreeNode *p, double tempq[2][16], double peng[2])
{
    int ii, i, j, k, indx;
    double sl[4], pt_comp[2][16];
    
    s_ComputeCoeffsCoulomb(p);
    for (indx = 0; indx < 16; indx++) {
        peng[0] = 0.0;
        ii = 0;
        for (i = 0; i < s_order+1; i++) {
            for (j = 0; j < s_order+1-i; j++) {
                for (k = 0; k < s_order+1-i-j; k++) {
                    peng[0] += s_der_coeff[i][j][k][indx]
                             * s_a[i+2+s_kk[indx][0]]
                                  [j+2+s_kk[indx][1]]
                                  [k+2+s_kk[indx][2]]
                             * p->ms[indx][ii];
                    ii++;
                }
            }
        }
    
        pt_comp[0][indx] = tempq[0][indx] * peng[0];
    }

    s_ComputeCoeffs(p);
    for (indx = 0; indx < 16; indx++) {
        peng[1] = 0.0;
        ii = 0;
        for (i = 0; i < s_order+1; i++) {
            for (j = 0; j < s_order+1-i; j++) {
                for (k = 0; k < s_order+1-i-j; k++) {
                    peng[1] += s_der_coeff[i][j][k][indx]
                             * s_a[i+2+s_kk[indx][0]]
                                  [j+2+s_kk[indx][1]]
                                  [k+2+s_kk[indx][2]]
                             * p->ms[indx][ii];
                    ii++;
                }
            }
        }
        
        pt_comp[1][indx] = tempq[1][indx] * peng[1];
    }

    sl[0] = pt_comp[0][0] - pt_comp[1][0];
    sl[1] = s_eps * (pt_comp[1][1] + pt_comp[1][2] + pt_comp[1][3])
                  - (pt_comp[0][1] + pt_comp[0][2] + pt_comp[0][3]);
    sl[2] = - (pt_comp[0][4] + pt_comp[0][5] + pt_comp[0][6])
            + (pt_comp[1][4] + pt_comp[1][5] + pt_comp[1][6]) / s_eps;
    sl[3] = 0.0;

    for (i = 7; i < 16; i++) {
        sl[3] += pt_comp[1][i] - pt_comp[0][i];
    }

    peng[0] = sl[0] + sl[1];
    peng[1] = sl[2] + sl[3];

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeCoeffs(TreeNode *p)
{
/* COMP_TCOEFF computes the Taylor coefficients of the potential
 * using a recurrence formula. The center of the expansion is the
 * midpoint of the node P. s_target_position and s_order are globally defined. */
 
    double dx, dy, dz, ddx, ddy, ddz, dist, fac;
    double kappax, kappay, kappaz;
    int i, j, k;

  /* setup variables */

    dx = s_target_position[0] - p->x_mid;
    dy = s_target_position[1] - p->y_mid;
    dz = s_target_position[2] - p->z_mid;

    ddx = 2.0 * dx;
    ddy = 2.0 * dy;
    ddz = 2.0 * dz;

    kappax = s_kappa * dx;
    kappay = s_kappa * dy;
    kappaz = s_kappa * dz;

    dist = dx*dx + dy*dy + dz*dz;
    fac = 1.0/dist;
    dist = sqrt(dist);

  /* 0th coeff or function val */
    s_b[2][2][2] = exp(-s_kappa * dist);
    s_a[2][2][2] = s_b[2][2][2] / dist;

  /* 2 indices are 0 */

    s_b[3][2][2] = kappax * s_a[2][2][2];
    s_b[2][3][2] = kappay * s_a[2][2][2];
    s_b[2][2][3] = kappaz * s_a[2][2][2];

    s_a[3][2][2] = fac * dx * (s_a[2][2][2] + s_kappa * s_b[2][2][2]);
    s_a[2][3][2] = fac * dy * (s_a[2][2][2] + s_kappa * s_b[2][2][2]);
    s_a[2][2][3] = fac * dz * (s_a[2][2][2] + s_kappa * s_b[2][2][2]);

    for (i = 2; i < s_order+3; i++) {
        s_b[i+2][2][2] = s_cf1[i-1] * s_kappa * (dx * s_a[i+1][2][2] - s_a[i][2][2]);
        s_b[2][i+2][2] = s_cf1[i-1] * s_kappa * (dy * s_a[2][i+1][2] - s_a[2][i][2]);
        s_b[2][2][i+2] = s_cf1[i-1] * s_kappa * (dz * s_a[2][2][i+1] - s_a[2][2][i]);

        s_a[i+2][2][2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][2][2]
                                  - s_cf3[i-1] * s_a[i][2][2]
                     + s_cf1[i-1] * s_kappa * (dx * s_b[i+1][2][2]
                                              - s_b[i][2][2]));

        s_a[2][i+2][2] = fac * (ddy * s_cf2[i-1] * s_a[2][i+1][2]
                                  - s_cf3[i-1] * s_a[2][i][2]
                     + s_cf1[i-1] * s_kappa * (dy * s_b[2][i+1][2]
                                              - s_b[2][i][2]));
        
        s_a[2][2][i+2] = fac * (ddz * s_cf2[i-1] * s_a[2][2][i+1]
                                  - s_cf3[i-1] * s_a[2][2][i]
                     + s_cf1[i-1] * s_kappa * (dz * s_b[2][2][i+1]
                                              - s_b[2][2][i]));
    }

  /* 1 index 0, 1 index 1, other >=1 */
    s_b[3][3][2] = kappax * s_a[2][3][2];
    s_b[3][2][3] = kappax * s_a[2][2][3];
    s_b[2][3][3] = kappay * s_a[2][2][3];

    s_a[3][3][2] = fac * (dx * s_a[2][3][2] + ddy * s_a[3][2][2] + kappax * s_b[2][3][2]);
    s_a[3][2][3] = fac * (dx * s_a[2][2][3] + ddz * s_a[3][2][2] + kappax * s_b[2][2][3]);
    s_a[2][3][3] = fac * (dy * s_a[2][2][3] + ddz * s_a[2][3][2] + kappay * s_b[2][2][3]);

    for (i = 2; i < s_order+2; i++) {
        s_b[3][2][i+2] = kappax * s_a[2][2][i+2];
        s_b[2][3][i+2] = kappay * s_a[2][2][i+2];
        s_b[2][i+2][3] = kappaz * s_a[2][i+2][2];
        s_b[3][i+2][2] = kappax * s_a[2][i+2][2];
        s_b[i+2][3][2] = kappay * s_a[i+2][2][2];
        s_b[i+2][2][3] = kappaz * s_a[i+2][2][2];

        s_a[3][2][i+2]=fac * (dx * s_a[2][2][i+2] + ddz * s_a[3][2][i+1] - s_a[3][2][i]
                      + kappax * s_b[2][2][i+2]);
        s_a[2][3][i+2]=fac * (dy * s_a[2][2][i+2] + ddz * s_a[2][3][i+1] - s_a[2][3][i]
                      + kappay * s_b[2][2][i+2]);
        s_a[2][i+2][3]=fac * (dz * s_a[2][i+2][2] + ddy * s_a[2][i+1][3] - s_a[2][i][3]
                      + kappaz * s_b[2][i+2][2]);
        s_a[3][i+2][2]=fac * (dx * s_a[2][i+2][2] + ddy * s_a[3][i+1][2] - s_a[3][i][2]
                      + kappax * s_b[2][i+2][2]);
        s_a[i+2][3][2]=fac * (dy * s_a[i+2][2][2] + ddx * s_a[i+1][3][2] - s_a[i][3][2]
                      + kappay * s_b[i+2][2][2]);
        s_a[i+2][2][3]=fac * (dz * s_a[i+2][2][2] + ddx * s_a[i+1][2][3] - s_a[i][2][3]
                      + kappaz * s_b[i+2][2][2]);
    }

    /* 1 index 0, others >=2 */
    for (i = 2; i < s_order+1; i++) {
        for (j = 2; j < s_order+3-i; j++) {
            s_b[i+2][j+2][2] = s_cf1[i-1] * s_kappa * (dx * s_a[i+1][j+2][2] - s_a[i][j+2][2]);
            s_b[i+2][2][j+2] = s_cf1[i-1] * s_kappa * (dx * s_a[i+1][2][j+2] - s_a[i][2][j+2]);
            s_b[2][i+2][j+2] = s_cf1[i-1] * s_kappa * (dy * s_a[2][i+1][j+2] - s_a[2][i][j+2]);

            s_a[i+2][j+2][2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][j+2][2]
                                             + ddy * s_a[i+2][j+1][2]
                                        - s_cf3[i-1] * s_a[i][j+2][2]
                                                   - s_a[i+2][j][2]
                           + s_cf1[i-1] * s_kappa * (dx * s_b[i+1][j+2][2]
                                                    - s_b[i][j+2][2]));

            s_a[i+2][2][j+2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][2][j+2]
                                             + ddz * s_a[i+2][2][j+1]
                                        - s_cf3[i-1] * s_a[i][2][j+2]
                                                   - s_a[i+2][2][j]
                           + s_cf1[i-1] * s_kappa * (dx * s_b[i+1][2][j+2]
                                                    - s_b[i][2][j+2]));

            s_a[2][i+2][j+2] = fac * (ddy * s_cf2[i-1] * s_a[2][i+1][j+2]
                                             + ddz * s_a[2][i+2][j+1]
                                        - s_cf3[i-1] * s_a[2][i][j+2]
                                                   - s_a[2][i+2][j]
                           + s_cf1[i-1] * s_kappa * (dy * s_b[2][i+1][j+2]
                                                    - s_b[2][i][j+2]));
        }
    }

  /* 2 indices 1,other >= 1 */
    s_b[3][3][3] = kappax * s_a[2][3][3];
    s_a[3][3][3] = fac * (dx * s_a[2][3][3]
                     + ddy * s_a[3][2][3]
                     + ddz * s_a[3][3][2]
                  + kappax * s_b[2][3][3]);

    for (i = 2; i < s_order+1; i++) {
        s_b[3][3][i+2] = kappax * s_a[2][3][i+2];
        s_b[3][i+2][3] = kappax * s_a[2][i+2][3];
        s_b[i+2][3][3] = kappay * s_a[i+2][2][3];

        s_a[3][3][i+2] = fac * (dx * s_a[2][3][i+2]
                           + ddy * s_a[3][2][i+2]
                           + ddz * s_a[3][3][i+1]
                                  -s_a[3][3][i]
                        + kappax * s_b[2][3][i+2]);

        s_a[3][i+2][3] = fac * (dx * s_a[2][i+2][3]
                           + ddy * s_a[3][i+1][3]
                           + ddz * s_a[3][i+2][2]
                                 - s_a[3][i][3]
                        + kappax * s_b[2][i+2][3]);

        s_a[i+2][3][3] = fac * (dy * s_a[i+2][2][3]
                           + ddx * s_a[i+1][3][3]
                           + ddz * s_a[i+2][3][2]
                                 - s_a[i][3][3]
                        + kappay * s_b[i+2][2][3]);
    }

  /* 1 index 1, others >=2 */
    for (i = 2; i < s_order; i++) {
        for (j = 2; j < s_order+2-i; j++) {
            s_b[3][i+2][j+2] = kappax * s_a[2][i+2][j+2];
            s_b[i+2][3][j+2] = kappay * s_a[i+2][2][j+2];
            s_b[i+2][j+2][3] = kappaz * s_a[i+2][j+2][2];

            s_a[3][i+2][j+2] = fac * (dx * s_a[2][i+2][j+2]
                                 + ddy * s_a[3][i+1][j+2]
                                 + ddz * s_a[3][i+2][j+1]
                                       - s_a[3][i][j+2]
                                       - s_a[3][i+2][j]
                              + kappax * s_b[2][i+2][j+2]);

            s_a[i+2][3][j+2] = fac * (dy * s_a[i+2][2][j+2]
                                 + ddx * s_a[i+1][3][j+2]
                                 + ddz * s_a[i+2][3][j+1]
                                       - s_a[i][3][j+2]
                                       - s_a[i+2][3][j]
                              + kappay * s_b[i+2][2][j+2]);
            
            s_a[i+2][j+2][3] = fac * (dz * s_a[i+2][j+2][2]
                                 + ddx * s_a[i+1][j+2][3]
                                 + ddy * s_a[i+2][j+1][3]
                                       - s_a[i][j+2][3]
                                       - s_a[i+2][j][3]
                              + kappaz * s_b[i+2][j+2][2]);
        }
    }

  /* all indices >=2 */
    for (k = 2; k < s_order-1; k++) {
        for (j = 2; j < s_order+1-k; j++) {
            for (i = 2; i < s_order+3-k-j; i++) {
                /*
                s_b[i+2][j+2][k+2] = s_cf1[i-1] * s_kappa * (dx * s_a[i+1][j+2][k+2]
                                                          - s_a[i][j+2][k+2]);

                s_a[i+2][j+2][k+2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][j+2][k+2]
                                                   + ddy * s_a[i+2][j+1][k+2]
                                                   + ddz * s_a[i+2][j+2][k+1]
                                              - s_cf3[i-1] * s_a[i][j+2][k+2]
                                                         - s_a[i+2][j][k+2]
                                                         - s_a[i+2][j+2][k]
                                 + s_cf1[i-1] * s_kappa * (dx * s_b[i+1][j+2][k+2]
                                                          - s_b[i][j+2][k+2]));
                */

                s_b[i+2][j+2][k+2] = s_cf1[i+j+k-1] * s_kappa 
                                                    * (dx * s_a[i+1][j+2][k+2]
                                                     + dy * s_a[i+2][j+1][k+2] 
                                                     + dz * s_a[i+2][j+2][k+1] 
                                                          - s_a[i][j+2][k+2]
                                                          - s_a[i+2][j][k+2]
                                                          - s_a[i+2][j+2][k]);

                s_a[i+2][j+2][k+2] = fac
                                   * (s_cf2[i+j+k-1] * (ddx * s_a[i+1][j+2][k+2]
                                                      + ddy * s_a[i+2][j+1][k+2]
                                                      + ddz * s_a[i+2][j+2][k+1])
                                    - s_cf3[i+j+k-1] * (s_a[i][j+2][k+2]
                                                      + s_a[i+2][j][k+2]
                                                      + s_a[i+2][j+2][k])
                                    + s_cf1[i+j+k-1] * s_kappa 
                                                     * (dx * s_b[i+1][j+2][k+2]
                                                     +  dy * s_b[i+2][j+1][k+2]
                                                     +  dz * s_b[i+2][j+2][k+1]
                                                           - s_b[i][j+2][k+2]
                                                           - s_b[i+2][j][k+2]
                                                           - s_b[i+2][j+2][k]));
            }
        }
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeCoeffsCoulomb(TreeNode *p)
{
/* COMP_TCOEFF computes the Taylor coefficients of the potential
 * using a recurrence formula. The center of the expansion is the
 * midpoint of the node P. s_target_position and s_order are globally defined. */
    double dx, dy, dz, ddx, ddy, ddz, dist, fac;
    int i, j, k;

  /* setup variables */

    dx = s_target_position[0] - p->x_mid;
    dy = s_target_position[1] - p->y_mid;
    dz = s_target_position[2] - p->z_mid;

    ddx = 2.0 * dx;
    ddy = 2.0 * dy;
    ddz = 2.0 * dz;

    dist = dx*dx + dy*dy + dz*dz;
    fac = 1.0 / dist;
    dist = sqrt(dist);

  /* 0th coeff or function val */
    s_a[2][2][2] = 1.0 / dist;

  /* 2 indices are 0 */

    s_a[3][2][2] = fac * dx * s_a[2][2][2];
    s_a[2][3][2] = fac * dy * s_a[2][2][2];
    s_a[2][2][3] = fac * dz * s_a[2][2][2];

    for (i = 2; i < s_order+3; i++) {

        s_a[i+2][2][2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][2][2]
                                    - s_cf3[i-1] * s_a[i][2][2]);

        s_a[2][i+2][2] = fac * (ddy * s_cf2[i-1] * s_a[2][i+1][2]
                                    - s_cf3[i-1] * s_a[2][i][2]);
        
        s_a[2][2][i+2] = fac * (ddz * s_cf2[i-1] * s_a[2][2][i+1]
                                    - s_cf3[i-1] * s_a[2][2][i]);
    }

  /* 1 index 0, 1 index 1, other >=1 */
    s_a[3][3][2] = fac * (dx * s_a[2][3][2] + ddy * s_a[3][2][2]);
    s_a[3][2][3] = fac * (dx * s_a[2][2][3] + ddz * s_a[3][2][2]);
    s_a[2][3][3] = fac * (dy * s_a[2][2][3] + ddz * s_a[2][3][2]);

    for (i = 2; i < s_order+2; i++) {

        s_a[3][2][i+2]=fac * (dx * s_a[2][2][i+2] + ddz * s_a[3][2][i+1] - s_a[3][2][i]);
        s_a[2][3][i+2]=fac * (dy * s_a[2][2][i+2] + ddz * s_a[2][3][i+1] - s_a[2][3][i]);
        s_a[2][i+2][3]=fac * (dz * s_a[2][i+2][2] + ddy * s_a[2][i+1][3] - s_a[2][i][3]);
        s_a[3][i+2][2]=fac * (dx * s_a[2][i+2][2] + ddy * s_a[3][i+1][2] - s_a[3][i][2]);
        s_a[i+2][3][2]=fac * (dy * s_a[i+2][2][2] + ddx * s_a[i+1][3][2] - s_a[i][3][2]);
        s_a[i+2][2][3]=fac * (dz * s_a[i+2][2][2] + ddx * s_a[i+1][2][3] - s_a[i][2][3]);
    }

    /* 1 index 0, others >=2 */
    for (i = 2; i < s_order+1; i++) {
        for (j = 2; j < s_order+3-i; j++) {

            s_a[i+2][j+2][2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][j+2][2]
                                                 + ddy * s_a[i+2][j+1][2]
                                          - s_cf3[i-1] * s_a[i][j+2][2]
                                                       - s_a[i+2][j][2]);

            s_a[i+2][2][j+2] = fac * (ddx * s_cf2[i-1] * s_a[i+1][2][j+2]
                                                 + ddz * s_a[i+2][2][j+1]
                                          - s_cf3[i-1] * s_a[i][2][j+2]
                                                       - s_a[i+2][2][j]);

            s_a[2][i+2][j+2] = fac * (ddy * s_cf2[i-1] * s_a[2][i+1][j+2]
                                                 + ddz * s_a[2][i+2][j+1]
                                          - s_cf3[i-1] * s_a[2][i][j+2]
                                                       - s_a[2][i+2][j]);
        }
    }

  /* 2 indices 1,other >= 1 */
    s_a[3][3][3] = fac * (dx * s_a[2][3][3]
                       + ddy * s_a[3][2][3]
                       + ddz * s_a[3][3][2]);

    for (i = 2; i < s_order+1; i++) {

        s_a[3][3][i+2] = fac * (dx * s_a[2][3][i+2]
                           + ddy * s_a[3][2][i+2]
                           + ddz * s_a[3][3][i+1]
                                  -s_a[3][3][i]);

        s_a[3][i+2][3] = fac * (dx * s_a[2][i+2][3]
                             + ddy * s_a[3][i+1][3]
                             + ddz * s_a[3][i+2][2]
                                   - s_a[3][i][3]);

        s_a[i+2][3][3] = fac * (dy * s_a[i+2][2][3]
                             + ddx * s_a[i+1][3][3]
                             + ddz * s_a[i+2][3][2]
                                   - s_a[i][3][3]);
    }

  /* 1 index 1, others >=2 */
    for (i = 2; i < s_order; i++) {
        for (j = 2; j < s_order+2-i; j++) {

            s_a[3][i+2][j+2] = fac * (dx * s_a[2][i+2][j+2]
                                   + ddy * s_a[3][i+1][j+2]
                                   + ddz * s_a[3][i+2][j+1]
                                         - s_a[3][i][j+2]
                                         - s_a[3][i+2][j]);

            s_a[i+2][3][j+2] = fac * (dy * s_a[i+2][2][j+2]
                                   + ddx * s_a[i+1][3][j+2]
                                   + ddz * s_a[i+2][3][j+1]
                                         - s_a[i][3][j+2]
                                         - s_a[i+2][3][j]);
            
            s_a[i+2][j+2][3] = fac * (dz * s_a[i+2][j+2][2]
                                   + ddx * s_a[i+1][j+2][3]
                                   + ddy * s_a[i+2][j+1][3]
                                         - s_a[i][j+2][3]
                                         - s_a[i+2][j][3]);
        }
    }

  /* all indices >=2 */
    for (k = 2; k < s_order-1; k++) {
        for (j = 2; j < s_order+1-k; j++) {
            for (i = 2; i < s_order+3-k-j; i++) {



                /* 
                s_a[i+2][j+2][k+2] = fac
                                   * (ddx * s_cf2[i-1] * s_a[i+1][j+2][k+2]
                                                 + ddy * s_a[i+2][j+1][k+2]
                                                 + ddz * s_a[i+2][j+2][k+1]
                                          - s_cf3[i-1] * s_a[i][j+2][k+2]
                                                       - s_a[i+2][j][k+2]
                                                       - s_a[i+2][j+2][k]);
                */

                s_a[i+2][j+2][k+2] = fac
                                   * (s_cf2[i+j+k-1] * (ddx * s_a[i+1][j+2][k+2]
                                                      + ddy * s_a[i+2][j+1][k+2]
                                                      + ddz * s_a[i+2][j+2][k+1])
                                    - s_cf3[i+j+k-1] * (s_a[i][j+2][k+2]
                                                      + s_a[i+2][j][k+2]
                                                      + s_a[i+2][j+2][k]));
            }
        }
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeDirectPB(int ibeg, int iend,
                             double *tpoten_old, double peng[2])
{
  /* COMPF_DIRECT directly computes the force on the current target
 * particle determined by the global variable s_target_position.*/
    int j;
    double peng_old[2], L1, L2, L3, L4, area;
    double tp[3], tq[3], sp[3], sq[3], r_s[3];
    double rs, irs, sumrs;
    double G0, kappa_rs, exp_kappa_rs, Gk;
    double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
    double G10, G20, G1, G2, G3, G4;

    peng[0] = 0.0;
    peng[1] = 0.0;

    tp[0] = s_target_position[0];
    tp[1] = s_target_position[1];
    tp[2] = s_target_position[2];
    tq[0] = s_target_normal[0];
    tq[1] = s_target_normal[1];
    tq[2] = s_target_normal[2];
    
    for (j = ibeg; j < iend+1; j++) {
        sp[0] = s_particle_position[0][j];
        sp[1] = s_particle_position[1][j];
        sp[2] = s_particle_position[2][j];
        sq[0] = s_particle_normal[0][j];
        sq[1] = s_particle_normal[1][j];
        sq[2] = s_particle_normal[2][j];

        r_s[0] = sp[0]-tp[0];  r_s[1] = sp[1]-tp[1];  r_s[2] = sp[2]-tp[2];
        sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
        rs = sqrt(sumrs);
        irs = 1.0/rs;
        G0 = ONE_OVER_4PI * irs;
        kappa_rs = s_kappa * rs;
        exp_kappa_rs = exp(-kappa_rs);
        Gk = exp_kappa_rs * G0;

        cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
        cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
        tp1 = G0* irs;
        tp2 = (1.0 + kappa_rs) * exp_kappa_rs;

        G10 = cos_theta0 * tp1;
        G20 = tp2 * G10;

        G1 = cos_theta * tp1;
        G2 = tp2 * G1;

        dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
        G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
        G4 = tp2*G3 - s_kappa2*cos_theta0*cos_theta*Gk;

        L1 = G1 - s_eps*G2;
        L2 = G0 - Gk;
        L3 = G4 - G3;
        L4 = G10 - G20/s_eps;

        peng_old[0] = tpoten_old[j];  
        peng_old[1] = tpoten_old[j + s_numpars];
        area = s_particle_area[j];

        peng[0] += (L1*peng_old[0] + L2*peng_old[1]) * area;
        peng[1] += (L3*peng_old[0] + L4*peng_old[1]) * area;
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_RemoveMoments(TreeNode *p)
{
/* REMOVE_NODE recursively removes each node from the
 * tree and deallocates its memory for MS array if it exits. */
 
    int i;

    if (p->exist_ms == 1) {
        free_matrix(p->ms);
        p->exist_ms = 0;
    }

    if (p->num_children > 0) {
        for (i = 0; i < p->num_children; i++)
            s_RemoveMoments(p->child[i]);
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_RemoveNode(TreeNode *p)
{
/* REMOVE_NODE recursively removes each node from the
 * tree and deallocates its memory for MS array if it exits. */
 
    int i;

    if (p->num_children > 0) {
        for (i = 0; i < 8; i++) {
            s_RemoveNode(p->child[i]);
            free(p->child[i]);
        }
        free(p->child);
    }

    return 0;
}
/**********************************************************/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* PRIVATE functions for preconditioning                     * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
static void leaflength(TreeNode *p, int idx, int *nrow)
{
/* find the leaf length */

    int i;

    if (idx == p->ibeg && p->num_children == 0) {
        *nrow = p->numpar;
    } else {
        if (p->num_children != 0) {
            for (i = 0; i < p->num_children; i++)
                leaflength(p->child[i], idx, nrow);
        }
     }

}
/**********************************************************/


/**********************************************************/
static int lu_decomp(double **A, int N, int *ipiv)
{
/* we're doing it this way because something is wrong with
 * linking with CMake */
 
  int i, j, k, imax;
  double maxA, *ptr, absA, Tol = 1.0e-14;

  for ( i = 0; i <= N; i++ )
    ipiv[i] = i; // record pivoting number

  for ( i = 0; i < N; i++ ) {
    maxA = 0.0;
    imax = i;
    for (k = i; k < N; k++)
      if ((absA = fabs(A[k][i])) > maxA) {
        maxA = absA;
        imax = k;
      }

    if (maxA < Tol) return 0; //failure, matrix is degenerate

    if (imax != i) {
      //pivoting P
      j = ipiv[i];
      ipiv[i] = ipiv[imax];
      ipiv[imax] = j;

      //pivoting rows of A
      ptr = A[i];
      A[i] = A[imax];
      A[imax] = ptr;

      //counting pivots starting from N (for determinant)
      ipiv[N]++;
    }

    for (j = i + 1; j < N; j++) {
      A[j][i] /= A[i][i];

      for (k = i + 1; k < N; k++)
        A[j][k] -= A[j][i] * A[i][k];
    }
  }

  return 1;
}
/**********************************************************/


/**********************************************************/
static void lu_solve(double **matrixA, int N, int *ipiv, double *rhs)
{
  /* b will contain the solution */
  
  int i, k;
  double *xtemp;

  make_vector(xtemp, N);

  for (i = 0; i < N; i++) {
    xtemp[i] = rhs[ipiv[i]];

    for (k = 0; k < i; k++)
      xtemp[i] -= matrixA[i][k] * xtemp[k];
  }

  for (i = N - 1; i >= 0; i--) {
    for (k = i + 1; k < N; k++)
      xtemp[i] -= matrixA[i][k] * xtemp[k];

    xtemp[i] = xtemp[i] / matrixA[i][i];
  }

  for (i = 0; i < N; i++) {
    rhs[i] = xtemp[i];
  }
  free_vector(xtemp);
}
/**********************************************************/


