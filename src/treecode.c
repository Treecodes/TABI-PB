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
#include <float.h>

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

#include "tree/struct_tree_linked_list_node.h"
#include "tree/tree_linked_list.h"

#include "struct_particles.h"
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
static int s_torder_lim;
static int s_torder3;

/* variable used by kernel independent moment computation */
double *tt, *ww;

/* these point to arrays located in TreeParticles */
static double *s_particle_position_x = NULL;
static double *s_particle_position_y = NULL;
static double *s_particle_position_z = NULL;

static double *s_particle_normal_x = NULL;
static double *s_particle_normal_y = NULL;
static double *s_particle_normal_z = NULL;

static double *s_particle_area = NULL;
static double *s_source_term = NULL;

/* global variables used when computing potential/force */
static double s_target_position[3];
static double s_target_normal[3];

static double **s_target_charge = NULL;
static double **s_source_charge = NULL;

/* global variables for reordering arrays */
//static int *s_order_arr = NULL;

/* root node of tree */
static struct TreeLinkedListNode *s_tree_root = NULL;


/* internal functions */
static int s_Setup(double xyz_limits[6], struct Particles *particles);

//static int s_CreateTree(struct TreeLinkedListNode *p, int ibeg, int iend, double xyzmm[6],
//                        int level);
//static int s_PartitionEight(double xyzmms[6][8], double xl, double yl,
//                            double zl, double lmax, double x_mid, double y_mid,
//                            double z_mid, int ind[8][2]);
                            
static int s_ComputePBKernel(double *phi);
static int s_ComputeAllMoments(struct TreeLinkedListNode *p, int ifirst);
static int s_ComputeMoments(struct TreeLinkedListNode *p);
static int s_RunTreecode(struct TreeLinkedListNode *p, double *tpoten_old,
                         double tempq[16], double peng[2]);
static int s_ComputeTreePB(struct TreeLinkedListNode *p, double tempq[16], double peng[2]);
static int s_ComputeDirectPB(int ibeg, int iend, double *tpoten_old,
                             double peng[2]);
                             
static int s_RemoveMoments(struct TreeLinkedListNode *p);
//static int s_RemoveNode(struct TreeLinkedListNode *p);

/* internal preconditioning functions */
static void leaflength(struct TreeLinkedListNode *p, int idx, int *nrow);
static int lu_decomp(double **A, int N, int *ipiv);
static void lu_solve(double **matrixA, int N, int *ipiv, double *rhs);



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* TreecodeInitialization and Finalization are used by       * * * */
/* tabipb() to interface with the treecode                   * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
int TreecodeInitialization(TABIPBparm *parm, struct Particles *particles)
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
    s_numpars = particles->num;
    s_order = parm->order;
    s_max_per_leaf = parm->maxparnode;
    theta = parm->theta;
    
    s_kappa = parm->kappa;
    s_kappa2 = parm->kappa2;
    s_eps = parm->eps;
    
    s_torder_lim = s_order+1;
    s_torder3 = s_torder_lim * s_torder_lim * s_torder_lim;
    
    s_min_level = 50000;
    s_max_level = 0;

    level = 0;

    
    s_particle_position_x = particles->x;
    s_particle_position_y = particles->y;
    s_particle_position_z = particles->z;
    
    s_particle_normal_x = particles->nx;
    s_particle_normal_y = particles->ny;
    s_particle_normal_z = particles->nz;
    
    s_particle_area = particles->area;
    s_source_term = particles->source_term;

    make_matrix(temp_normal, 3, s_numpars);
    make_vector(temp_area, s_numpars);
    make_vector(temp_source, 2 * s_numpars);
    

/* Call SETUP to allocate arrays for Taylor expansions */
/* and setup global variables. Also, copy variables into global copy arrays. */
    s_Setup(xyz_limits, particles);

    //s_tree_root = (TreeNode*)calloc(1, sizeof(TreeNode));
    //s_CreateTree(s_tree_root, 0, s_numpars-1, xyz_limits, level);
    
    int numnodes, numleaves, max_depth;
    
    TreeLinkedList_Construct(&s_tree_root, NULL, particles, 0, particles->num-1,
                s_max_per_leaf, xyz_limits, &numnodes, &numleaves,
                &s_min_level, &s_max_level, &max_depth, 0);
    
    if (rank == 0) {
        printf("Created tree for %d particles with max %d per node.\n\n",
               s_numpars, s_max_per_leaf);
    }

    memcpy(temp_normal[0], s_particle_normal_x, s_numpars*sizeof(double));
    memcpy(temp_normal[1], s_particle_normal_y, s_numpars*sizeof(double));
    memcpy(temp_normal[2], s_particle_normal_z, s_numpars*sizeof(double));
    memcpy(temp_area, s_particle_area, s_numpars*sizeof(double));
    memcpy(temp_source, s_source_term, 2*s_numpars*sizeof(double));
    
    for (i = 0; i < s_numpars; i++) {
        s_particle_normal_x[i]    = temp_normal[0][particles->order[i]];
        s_particle_normal_y[i]    = temp_normal[1][particles->order[i]];
        s_particle_normal_z[i]    = temp_normal[2][particles->order[i]];
        s_particle_area[i]           =   temp_area[particles->order[i]];
        s_source_term[i]             = temp_source[particles->order[i]];
        s_source_term[i + s_numpars] = temp_source[particles->order[i] + s_numpars];
    }

    free_matrix(temp_normal);
    free_vector(temp_area);
    free_vector(temp_source);

    make_matrix(s_target_charge, s_numpars, 16);
    make_matrix(s_source_charge, s_numpars, 16);

    return 0;
}
/**********************************************************/


/********************************************************/
int TreecodeFinalization(struct Particles *particles)
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

    memcpy(temp_position[0], particles->x, s_numpars*sizeof(double));
    memcpy(temp_position[1], particles->y, s_numpars*sizeof(double));
    memcpy(temp_position[2], particles->z, s_numpars*sizeof(double));
    memcpy(temp_normal[0], particles->nx, s_numpars*sizeof(double));
    memcpy(temp_normal[1], particles->ny, s_numpars*sizeof(double));
    memcpy(temp_normal[2], particles->nz, s_numpars*sizeof(double));
    memcpy(temp_area, particles->area, s_numpars*sizeof(double));
    memcpy(temp_source, particles->source_term, 2*s_numpars*sizeof(double));
    memcpy(temp_xvct, particles->xvct, 2*s_numpars*sizeof(double));
    
    for (i = 0; i < s_numpars; i++) {
        particles->x[particles->order[i]]     = temp_position[0][i];
        particles->y[particles->order[i]]     = temp_position[1][i];
        particles->z[particles->order[i]]     = temp_position[2][i];
        particles->nx[particles->order[i]]    = temp_normal[0][i];
        particles->ny[particles->order[i]]    = temp_normal[1][i];
        particles->nz[particles->order[i]]    = temp_normal[2][i];
        particles->area[particles->order[i]]         = temp_area[i];
        particles->source_term[particles->order[i]]  = temp_source[i];
        particles->source_term[particles->order[i] + s_numpars]
                                                     = temp_source[i + s_numpars];
        particles->xvct[particles->order[i]]         = temp_xvct[i];
        particles->xvct[particles->order[i] + s_numpars]
                                                     = temp_xvct[i + s_numpars];
    }

    free_matrix(temp_position);
    free_matrix(temp_normal);
    free_vector(temp_area);
    free_vector(temp_source);
    free_vector(temp_xvct);

/***********treecode_initialization*******/

    free_matrix(s_target_charge);
    free_matrix(s_source_charge);
    
/***********clean tree structure**********/

    TreeLinkedList_Free(&s_tree_root);

/***********variables in setup************/

    free_vector(particles->order);

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
    double temp_charge[16];
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
            s_target_position[0] = s_particle_position_x[i];
            s_target_position[1] = s_particle_position_y[i];
            s_target_position[2] = s_particle_position_z[i];
            s_target_normal[0] = s_particle_normal_x[i];
            s_target_normal[1] = s_particle_normal_y[i];
            s_target_normal[2] = s_particle_normal_z[i];
        
            for (k = 0; k < 16; k++) {
                temp_charge[k] = s_target_charge[i][k];
            }

      /* remove the singularity */
            temp_x = s_particle_position_x[i];
            temp_area = s_particle_area[i];
            s_particle_position_x[i] += 100.123456789;
            s_particle_area[i] = 0.0;

      /* start to use Treecode */
            s_RunTreecode(s_tree_root, tpoten_old, temp_charge, peng);

            tpoten[i] = tpoten_temp[i] * *beta
                      + (pre1 * peng_old[0] - peng[0]) * *alpha;
            tpoten[s_numpars+i] = tpoten_temp[s_numpars+i] * *beta
                                + (pre2 * peng_old[1] - peng[1]) * *alpha;

            s_particle_position_x[i] = temp_x;
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
            tp[0] = s_particle_position_x[i];
            tp[1] = s_particle_position_y[i];
            tp[2] = s_particle_position_z[i];
            tq[0] = s_particle_normal_x[i];
            tq[1] = s_particle_normal_y[i];
            tq[2] = s_particle_normal_z[i];

            for (j = ibeg; j < i; j++) {
                sp[0] = s_particle_position_x[j];
                sp[1] = s_particle_position_y[j];
                sp[2] = s_particle_position_z[j];
                sq[0] = s_particle_normal_x[j];
                sq[1] = s_particle_normal_y[j];
                sq[2] = s_particle_normal_z[j];

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
                sp[0] = s_particle_position_x[j];
                sp[1] = s_particle_position_y[j];
                sp[2] = s_particle_position_z[j];
                sq[0] = s_particle_normal_x[j];
                sq[1] = s_particle_normal_y[j];
                sq[2] = s_particle_normal_z[j];

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
static int s_Setup(double xyz_limits[6], struct Particles *particles)
{
/* SETUP allocates and initializes arrays needed for the Taylor expansion.
 Also, global variables are set and the Cartesian coordinates of
 the smallest box containing the particles is determined. The particle
 positions and charges are copied so that they can be restored upon exit.*/
 
    int ierr;
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif
    
    if (rank == 0) {
        printf("Setting up arrays for Taylor expansion...\n");
    }
    
    make_vector(tt, s_torder_lim);
    make_vector(ww, s_torder_lim);
    
    /* initializing array for Chev points */
    for (int i = 0; i < s_torder_lim; i++)
        tt[i] = cos(i * M_PI / s_order);
    
    ww[0] = 0.25 * (s_order*s_order/3.0 + 1.0/6.0);
    ww[s_order] = -ww[0];
    
    for (int i = 1; i < s_order; i++) {
        double xx = i * M_PI / s_order;
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
    }

/* find bounds of Cartesion box enclosing the particles */

    xyz_limits[0] = MinVal(s_particle_position_x, s_numpars);
    xyz_limits[1] = MaxVal(s_particle_position_x, s_numpars);
    xyz_limits[2] = MinVal(s_particle_position_y, s_numpars);
    xyz_limits[3] = MaxVal(s_particle_position_y, s_numpars);
    xyz_limits[4] = MinVal(s_particle_position_z, s_numpars);
    xyz_limits[5] = MaxVal(s_particle_position_z, s_numpars);

    //printf("x-limits of box: %f, %f\n", xyz_limits[0], xyz_limits[1]);
    //printf("y-limits of box: %f, %f\n", xyz_limits[2], xyz_limits[3]);
    //printf("z-limits of box: %f, %f\n", xyz_limits[4], xyz_limits[5]);

    make_vector(particles->order, particles->num);

    for (int i = 0; i < particles->num; i++) {
        particles->order[i] = i;
    }

    return 0;
}
/********************************************************/


/********************************************************/
//
//static int s_CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6],
//                        int level)
//{
///*CREATE_TREE recursively create the tree structure. Node P is
//  input, which contains particles indexed from IBEG to IEND. After
//  the node parameters are set subdivision occurs if IEND-IBEG+1 > s_max_per_leaf.
//  Real array XYZMM contains the min and max values of the coordinates
//  of the particle in P, thus defining the box. */
//
//  /* local variables */
//    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
//    int ind[8][2];
//    double xyzmms[6][8];
//    int i, j, loclev, numposchild;
//    double lxyzmm[6];
//
///* set node fields: number of particles, exist_ms and xyz bounds */
//    p->numpar = iend-ibeg+1;
//    p->exist_ms = 0;
//
//    p->x_min = xyzmm[0];
//    p->x_max = xyzmm[1];
//    p->y_min = xyzmm[2];
//    p->y_max = xyzmm[3];
//    p->z_min = xyzmm[4];
//    p->z_max = xyzmm[5];
//
///* compute aspect ratio */
//    xl = p->x_max-p->x_min;
//    yl = p->y_max-p->y_min;
//    zl = p->z_max-p->z_min;
//
//    lmax = xl;
//    if (lmax < yl) lmax = yl;
//    if (lmax < zl) lmax = zl;
//
//    t1 = lmax;
//    t2 = xl;
//    if (t2 > yl) t2 = yl;
//    if (t2 > zl) t2 = zl;
//
//    if (t2 != 0.0) {
//        p->aspect = t1/t2;
//    } else {
//        p->aspect = 0.0;
//    }
//
///* midpoint coordinates, RADIUS and SQRADIUS */
//    p->x_mid = (p->x_max + p->x_min) / 2.0;
//    p->y_mid = (p->y_max + p->y_min) / 2.0;
//    p->z_mid = (p->z_max + p->z_min) / 2.0;
//    t1 = p->x_max - p->x_mid;
//    t2 = p->y_max - p->y_mid;
//    t3 = p->z_max - p->z_mid;
//    p->radius = sqrt(t1*t1 + t2*t2 + t3*t3);
//
///* set particle limits, tree level of node, and nullify children pointers */
//    p->ibeg = ibeg;
//    p->iend = iend;
//    p->level = level;
//    if (s_max_level < level) s_max_level = level;
//
//    p->num_children = 0;
//
//    make_vector(p->child, 8);
//    for (i = 0; i < 8; i++) {
//        p->child[i] = (TreeNode*)calloc(1, sizeof(TreeNode));
//    }
//
//
//    if (p->numpar > s_max_per_leaf) {
///* set IND array to 0 and then call PARTITION routine. IND array holds indices
// * of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1 */
//
//        xyzmms[0][0] = p->x_min;
//        xyzmms[1][0] = p->x_max;
//        xyzmms[2][0] = p->y_min;
//        xyzmms[3][0] = p->y_max;
//        xyzmms[4][0] = p->z_min;
//        xyzmms[5][0] = p->z_max;
//
//        for (i = 0; i < 8; i++) {
//            ind[i][0] = 0;
//            ind[i][1] = 0;
//        }
//
//        ind[0][0] = ibeg;
//        ind[0][1] = iend;
//        x_mid = p->x_mid;
//        y_mid = p->y_mid;
//        z_mid = p->z_mid;
//
//        numposchild = s_PartitionEight(xyzmms, xl, yl, zl, lmax,
//                                       x_mid, y_mid, z_mid, ind);
//
///* Shrink the box */
//        for (i = 0; i < 8; i++) {
//            if (ind[i][0] < ind[i][1]) {
//                xyzmms[0][i] = MinVal(&s_particle_position_x[ind[i][0]],
//                                      ind[i][1]-ind[i][0]);
//                xyzmms[1][i] = MaxVal(&s_particle_position_x[ind[i][0]],
//                                      ind[i][1]-ind[i][0]);
//                xyzmms[2][i] = MinVal(&s_particle_position_y[ind[i][0]],
//                                      ind[i][1]-ind[i][0]);
//                xyzmms[3][i] = MaxVal(&s_particle_position_y[ind[i][0]],
//                                      ind[i][1]-ind[i][0]);
//                xyzmms[4][i] = MinVal(&s_particle_position_z[ind[i][0]],
//                                      ind[i][1]-ind[i][0]);
//                xyzmms[5][i] = MaxVal(&s_particle_position_z[ind[i][0]],
//                                      ind[i][1]-ind[i][0]);
//            }
//        }
///* create children if indicated and store info in parent */
//        loclev = level + 1;
//
//        for (i = 0; i < numposchild; i++) {
//            if (ind[i][0] <= ind[i][1]) {
//                p->num_children = p->num_children + 1;
//                for (j = 0; j < 6; j++) {
//                    lxyzmm[j] = xyzmms[j][i];
//                }
//                s_CreateTree(p->child[p->num_children-1],
//                             ind[i][0], ind[i][1], lxyzmm, loclev);
//            }
//        }
//    } else {
//        if (level < s_min_level) {
//            s_min_level = level;
//        }
//    }
//
//    return 0;
//}
///**********************************************************/
//
//
///********************************************************/
//static int s_PartitionEight(double xyzmms[6][8], double xl, double yl,
//                            double zl, double lmax, double x_mid, double y_mid,
//                            double z_mid, int ind[8][2])
//{
///* PARTITION_8 determines the particle indices of the eight sub boxes
// * containing the particles after the box defined by particles I_BEG
// * to I_END is divided by its midpoints in each coordinate direction.
// * The determination of the indices is accomplished by the subroutine
// * PARTITION. A box is divided in a coordinate direction as long as the
// * resulting aspect ratio is not too large. This avoids the creation of
// * "narrow" boxes in which Talyor expansions may become inefficient.
// * On exit the INTEGER array IND (dimension 8 x 2) contains
// * the indice limits of each new box (node) and NUMPOSCHILD the number
// * of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
// * that box J is empty.*/
//    int temp_ind, i, j;
//    double critlen;
//    int numposchild;
//
//    numposchild = 1;
//    critlen = lmax/sqrt(2.0);
//
//    if (xl >= critlen) {
//        temp_ind = Partition(s_particle_position_x, s_particle_position_y,
//                             s_particle_position_z, s_order_arr,
//                             ind[0][0], ind[0][1], x_mid);
//        ind[1][0] = temp_ind+1;
//        ind[1][1] = ind[0][1];
//        ind[0][1] = temp_ind;
//        for (i = 0; i < 6; i++) {
//            xyzmms[i][1] = xyzmms[i][0];
//        }
//        xyzmms[1][0] = x_mid;
//        xyzmms[0][1] = x_mid;
//        numposchild *= 2;
//    }
//
//    if (yl >= critlen) {
//        for (i = 0; i < numposchild; i++) {
//            temp_ind = Partition(s_particle_position_y, s_particle_position_x,
//                                 s_particle_position_z, s_order_arr,
//                                 ind[i][0], ind[i][1], y_mid);
//            ind[numposchild+i][0] = temp_ind+1;
//            ind[numposchild+i][1] = ind[i][1];
//            ind[i][1] = temp_ind;
//            for (j = 0; j < 6; j++) {
//                xyzmms[j][numposchild+i] = xyzmms[j][i];
//            }
//            xyzmms[3][i] = y_mid;
//            xyzmms[2][numposchild+i] = y_mid;
//        }
//        numposchild *= 2;
//    }
//
//    if (zl >= critlen) {
//        for (i = 0; i < numposchild; i++) {
//            temp_ind = Partition(s_particle_position_z, s_particle_position_x,
//                                 s_particle_position_y, s_order_arr,
//                                 ind[i][0], ind[i][1], z_mid);
//            ind[numposchild+i][0] = temp_ind+1;
//            ind[numposchild+i][1] = ind[i][1];
//            ind[i][1] = temp_ind;
//            for (j = 0; j < 6; j++) {
//                xyzmms[j][numposchild+i] = xyzmms[j][i];
//            }
//            xyzmms[5][i] = z_mid;
//            xyzmms[4][numposchild+i] = z_mid;
//        }
//        numposchild *= 2;
//    }
//
//    return (numposchild);
//}
/********************************************************/


/********************************************************/
static int s_ComputePBKernel(double *phi)
{

    int i, iknl, ixyz, jxyz, indx;

    for (i = 0; i < s_numpars; i++) {
        indx = 0;
        s_target_charge[i][indx] = 1.0;
        s_source_charge[i][indx] = s_particle_area[i] * phi[s_numpars+i];

        for (iknl = 0; iknl < 2; iknl++) {
//            for (ixyz = 0; ixyz < 3; ixyz++) {
//                indx += 1;
//
//                s_target_charge[i][indx] = 1.0 * (1-iknl) + s_particle_normal[ixyz][i] * iknl;
//
//                s_source_charge[i][indx] = (s_particle_normal[ixyz][i] * (1-iknl) + 1.0 * iknl)
//                                         * s_particle_area[i] * phi[iknl*s_numpars+i];
                        
                indx += 1;
                s_target_charge[i][indx] = 1.0 * (1-iknl) + s_particle_normal_x[i] * iknl;
                s_source_charge[i][indx] = (s_particle_normal_x[i] * (1-iknl) + 1.0 * iknl)
                                          * s_particle_area[i] * phi[iknl*s_numpars+i];
                                         
                indx += 1;
                s_target_charge[i][indx] = 1.0 * (1-iknl) + s_particle_normal_y[i] * iknl;
                s_source_charge[i][indx] = (s_particle_normal_y[i] * (1-iknl) + 1.0 * iknl)
                                          * s_particle_area[i] * phi[iknl*s_numpars+i];
                                         
                indx += 1;
                s_target_charge[i][indx] = 1.0 * (1-iknl) + s_particle_normal_z[i] * iknl;
                s_source_charge[i][indx] = (s_particle_normal_z[i] * (1-iknl) + 1.0 * iknl)
                                          * s_particle_area[i] * phi[iknl*s_numpars+i];
//            }
        }

        for (ixyz = 0; ixyz < 3; ixyz++) {
//            for (jxyz = 0; jxyz < 3; jxyz++) {
//                indx += 1;
//                s_target_charge[i][indx] =  s_particle_normal[jxyz][i];
//                s_source_charge[i][indx] = -s_particle_normal[ixyz][i] * s_particle_area[i] * phi[i];
                
                indx += 1;
                s_target_charge[i][indx] =  s_particle_normal_x[i];
                s_source_charge[i][indx] = -s_particle_normal_x[i] * s_particle_area[i] * phi[i];
                
                indx += 1;
                s_target_charge[i][indx] =  s_particle_normal_y[i];
                s_source_charge[i][indx] = -s_particle_normal_y[i] * s_particle_area[i] * phi[i];
                
                indx += 1;
                s_target_charge[i][indx] =  s_particle_normal_z[i];
                s_source_charge[i][indx] = -s_particle_normal_z[i] * s_particle_area[i] * phi[i];
//            }
        }
    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeAllMoments(struct TreeLinkedListNode *p, int ifirst)
{
/* REMOVE_NODE recursively removes each node from the tree and deallocates
 * its memory for MS array if it exits. */

    int i;
    
    if (p->exist_ms == 0 && ifirst == 0) {
        make_matrix(p->ms, 16, s_torder3);
        make_vector(p->tx, s_torder_lim);
        make_vector(p->ty, s_torder_lim);
        make_vector(p->tz, s_torder_lim);
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
static int s_ComputeMoments(struct TreeLinkedListNode *p)
{
/* COMP_MS computes the moments for node P needed in the Taylor
 * approximation */

    int i, j, k1, k2, k3, kk;
    double dx, dy, dz;
    double x0, x1, y0, y1, z0, z1;
    double sumA1, sumA2, sumA3;
    double temp11, temp12, temp21, temp22;
    double mom1, mom2, mom3, mom4, mom5, mom6, mom7, mom8;
    double xx, yy, zz;
    double *xibeg, *yibeg, *zibeg, *qq;
    
    double Dd, dj[s_torder_lim];
    double a1i[s_torder_lim], a2j[s_torder_lim], a3k[s_torder_lim];
    double w1i[s_torder_lim];
    double summ[16][s_torder3];
    int a1exactind, a2exactind, a3exactind;
    
    
    for (i = 0; i < 16; i++) {
        for (j = 0; j < s_torder3; j++) {
            p->ms[i][j] = 0.0;
        }
    }
    
    xibeg = &(s_particle_position_x[p->ibeg]);
    yibeg = &(s_particle_position_y[p->ibeg]);
    zibeg = &(s_particle_position_z[p->ibeg]);
    
    x0 = p->x_min;
    x1 = p->x_max;
    y0 = p->y_min;
    y1 = p->y_max;
    z0 = p->z_min;
    z1 = p->z_max;
    
    for (i = 0; i < 16; i++) {
        for (j = 0; j < s_torder3; j++) {
            summ[i][j] = 0.0;
        }
    }
    
    for (i = 0; i < s_torder_lim; i++) {
        p->tx[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        p->ty[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        p->tz[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    
    dj[0] = 0.5;
    dj[s_order] = 0.5;
    for (j = 1; j < s_order; j++)
        dj[j] = 1.0;
    
    for (j = 0; j < s_torder_lim; j++)
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];

    for (i = 0; i < p->numpar; i++) {
    
        sumA1 = 0.0;
        sumA2 = 0.0;
        sumA3 = 0.0;
        
        a1exactind = -1;
        a2exactind = -1;
        a3exactind = -1;
    
        xx = xibeg[i];
        yy = yibeg[i];
        zz = zibeg[i];
        qq = s_source_charge[p->ibeg+i];
        
        for (j = 0; j < s_torder_lim; j++) {
            a1i[j] = w1i[j] / (xx - p->tx[j]);
            a2j[j] = w1i[j] / (yy - p->ty[j]);
            a3k[j] = w1i[j] / (zz - p->tz[j]);

            sumA1 += a1i[j];
            sumA2 += a2j[j];
            sumA3 += a3k[j];
            
            if (fabs(xx - p->tx[j]) < DBL_MIN) a1exactind = j;
            if (fabs(yy - p->ty[j]) < DBL_MIN) a2exactind = j;
            if (fabs(zz - p->tz[j]) < DBL_MIN) a3exactind = j;
        }
        
        if (a1exactind > -1) {
            sumA1 = 1.0;
            for (j = 0; j < s_torder_lim; j++) a1i[j] = 0.0;
            a1i[a1exactind] = 1.0;
        }
        
        if (a2exactind > -1) {
            sumA2 = 1.0;
            for (j = 0; j < s_torder_lim; j++) a2j[j] = 0.0;
            a2j[a2exactind] = 1.0;
        }
        
        if (a3exactind > -1) {
            sumA3 = 1.0;
            for (j = 0; j < s_torder_lim; j++) a3k[j] = 0.0;
            a3k[a3exactind] = 1.0;
        }
        
        Dd = 1.0 / (sumA1 * sumA2 * sumA3);
        
        kk = -1;
        for (k1 = 0; k1 < s_torder_lim; k1++) {
            for (k2 = 0; k2 < s_torder_lim; k2++) {
                for (k3 = 0; k3 < s_torder_lim; k3++) {
                    kk++;
                
                    mom1 = a1i[k1] * a2j[k2] * a3k[k3] * Dd;
                    
                    for (j = 0; j < 7; j++) {
                        summ[j][kk] += mom1 * qq[j];
                    }
                    
                    for (j = 7; j < 14; j+=3) {
                        summ[j][kk] += mom1 * qq[j];
                    }
                }
            }
        }
    }
    
    for (j = 0; j < 7; j++)
        memcpy(p->ms[j], summ[j], s_torder3*sizeof(double));
    
    for (j = 7; j < 10; j++)
        memcpy(p->ms[j], summ[7], s_torder3*sizeof(double));
    
    for (j = 10; j < 13; j++)
        memcpy(p->ms[j], summ[10], s_torder3*sizeof(double));
    
    for (j = 13; j < 16; j++)
        memcpy(p->ms[j], summ[13], s_torder3*sizeof(double));


    return 0;
}
/********************************************************/


/********************************************************/
static int s_RunTreecode(struct TreeLinkedListNode *p, double *tpoten_old, double tempq[16],
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
static int s_ComputeTreePB(struct TreeLinkedListNode *p, double tempq[16], double peng[2])
{
    double sl[4], pt_comp[16];
    
    for (int indx = 0; indx < 16; indx++) {
        pt_comp[indx] = 0;
    }

    int ii = 0;
    for (int i = 0; i < s_torder_lim; i++) {
        for (int j = 0; j < s_torder_lim; j++) {
            for (int k = 0; k < s_torder_lim; k++) {

                double dx = s_target_position[0] - p->tx[i];
                double dy = s_target_position[1] - p->ty[j];
                double dz = s_target_position[2] - p->tz[k];

                double r2    = dx*dx + dy*dy + dz*dz;
                double r     = sqrt(r2);
                double rinv  = 1.0 / r;
                double r3inv = rinv  * rinv * rinv;
                double r5inv = r3inv * rinv * rinv;

                double expkr   =  exp(-s_kappa * r);
                double d1term  =  r3inv * expkr * (1. + (s_kappa * r)); 
                double d1term1 = -r3inv + d1term * s_eps;
                double d1term2 = -r3inv + d1term / s_eps;
                double d2term  =  r5inv * (-3. + expkr * (3. + (3. * s_kappa * r) + (s_kappa * s_kappa * r2)));
                double d3term  =  r3inv * ( 1. - expkr * (1. + s_kappa * r));


                pt_comp[0]  += (ONE_OVER_4PI * tempq[0]  * p->ms[0][ii]  *  rinv * (1. - expkr));

                pt_comp[1]  += (ONE_OVER_4PI * tempq[1]  * p->ms[1][ii]  *  dx * d1term1);
                pt_comp[2]  += (ONE_OVER_4PI * tempq[2]  * p->ms[2][ii]  *  dy * d1term1);
                pt_comp[3]  += (ONE_OVER_4PI * tempq[3]  * p->ms[3][ii]  *  dz * d1term1);

                pt_comp[4]  += (ONE_OVER_4PI * tempq[4]  * p->ms[4][ii]  *  dx * d1term2);
                pt_comp[5]  += (ONE_OVER_4PI * tempq[5]  * p->ms[5][ii]  *  dy * d1term2);
                pt_comp[6]  += (ONE_OVER_4PI * tempq[6]  * p->ms[6][ii]  *  dz * d1term2);

                pt_comp[7]  += (ONE_OVER_4PI * tempq[7]  * p->ms[7][ii]  * (dx * dx * d2term + d3term));
                pt_comp[8]  += (ONE_OVER_4PI * tempq[8]  * p->ms[8][ii]  *  dx * dy * d2term);
                pt_comp[9]  += (ONE_OVER_4PI * tempq[9]  * p->ms[9][ii]  *  dx * dz * d2term);

                pt_comp[10] += (ONE_OVER_4PI * tempq[10] * p->ms[10][ii] *  dx * dy * d2term);
                pt_comp[11] += (ONE_OVER_4PI * tempq[11] * p->ms[11][ii] * (dy * dy * d2term + d3term));
                pt_comp[12] += (ONE_OVER_4PI * tempq[12] * p->ms[12][ii] *  dy * dz * d2term);

                pt_comp[13] += (ONE_OVER_4PI * tempq[13] * p->ms[13][ii] *  dx * dz * d2term);
                pt_comp[14] += (ONE_OVER_4PI * tempq[14] * p->ms[14][ii] *  dy * dz * d2term);
                pt_comp[15] += (ONE_OVER_4PI * tempq[15] * p->ms[15][ii] * (dz * dz * d2term + d3term));

                ii++;
            }
        }
    }
    

    sl[0] = pt_comp[0];
    sl[1] = pt_comp[1] + pt_comp[2] + pt_comp[3];
    sl[2] = pt_comp[4] + pt_comp[5] + pt_comp[6];

    sl[3] = 0.0;
    for (int i = 7; i < 16; i++) {
        sl[3] += pt_comp[i];
    }

    peng[0] = sl[0] + sl[1];
    peng[1] = sl[2] + sl[3];

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
        sp[0] = s_particle_position_x[j];
        sp[1] = s_particle_position_y[j];
        sp[2] = s_particle_position_z[j];
        sq[0] = s_particle_normal_x[j];
        sq[1] = s_particle_normal_y[j];
        sq[2] = s_particle_normal_z[j];

        r_s[0] = sp[0]-tp[0];  
        r_s[1] = sp[1]-tp[1];  
        r_s[2] = sp[2]-tp[2];
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
static int s_RemoveMoments(struct TreeLinkedListNode *p)
{
/* REMOVE_NODE recursively removes each node from the
 * tree and deallocates its memory for MS array if it exits. */
 
    int i;

    if (p->exist_ms == 1) {
        free_matrix(p->ms);
        free_vector(p->tx);
        free_vector(p->ty);
        free_vector(p->tz);
        p->exist_ms = 0;
    }

    if (p->num_children > 0) {
        for (i = 0; i < p->num_children; i++)
            s_RemoveMoments(p->child[i]);
    }

    return 0;
}
/********************************************************/


///********************************************************/
//static int s_RemoveNode(TreeNode *p)
//{
///* REMOVE_NODE recursively removes each node from the
// * tree and deallocates its memory for MS array if it exits. */
//
//    int i;
//
//    if (p->num_children > 0) {
//        for (i = 0; i < 8; i++) {
//            s_RemoveNode(p->child[i]);
//            free(p->child[i]);
//        }
//        free(p->child);
//    }
//
//    return 0;
//}
///**********************************************************/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* PRIVATE functions for preconditioning                     * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
static void leaflength(struct TreeLinkedListNode *p, int idx, int *nrow)
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


