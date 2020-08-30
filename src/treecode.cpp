/**************************************************************************
* FILE NAME: treecode.c
*
* PURPOSE: Contains all treecode-related functions and variables,
*          including treecode initialization and finalization functions
*          that interface with tabipb.c, and matrix-vector multiplication
*          and solve functions that interface with run_gmres.c
*
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI
*          Jiahui Chen, Southern Methodist University, Dallas, TX
*
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:
*          Weihua Geng, Southern Methodist University, Dallas, TX
*          Robery Krasny, University of Michigan, Ann Arbor, MI
*
**************************************************************************/

#include <vector>
#include <numeric>
#include <iostream>

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

static double *s_target_charge_   = NULL;
static double *s_target_charge_dx = NULL;
static double *s_target_charge_dy = NULL;
static double *s_target_charge_dz = NULL;
static double *s_source_charge_   = NULL;
static double *s_source_charge_dx = NULL;
static double *s_source_charge_dy = NULL;
static double *s_source_charge_dz = NULL;

/* global variables for reordering arrays */
//static int *s_order_arr = NULL;

/* root node of tree */
static struct TreeLinkedListNode *s_tree_root = NULL;


/* internal functions */
static int s_Setup(double xyz_limits[6], struct Particles *particles);
static int s_ComputePBKernel(double *phi);
static int s_ComputeAllMoments(struct TreeLinkedListNode *p, int ifirst);
static int s_ComputeMoments(struct TreeLinkedListNode *p);
static int s_RunTreecode(struct TreeLinkedListNode *p, double *tpoten_old,
                         double tempq[4], double peng[2]);
static int s_ComputeTreePB(struct TreeLinkedListNode *p, double tempq[4], double peng[2]);
static int s_ComputeDirectPB(int ibeg, int iend, double *tpoten_old,
                             double peng[2]);
                             
static int s_RemoveMoments(struct TreeLinkedListNode *p);

/* internal preconditioning functions */
static void leaflength(struct TreeLinkedListNode *p, int idx, int *nrow);


template <typename order_iterator, typename value_iterator>
void reorder(order_iterator order_begin, order_iterator order_end, value_iterator v_begin)
{
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    
    auto v_end = v_begin + std::distance(order_begin, order_end) + 1;
    std::vector<value_t> tmp(v_begin, v_end);

    std::for_each(order_begin, order_end,
                  [&tmp, &v_begin](index_t idx){ *v_begin = tmp[idx]; std::advance(v_begin, 1); });
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* TreecodeInitialization and Finalization are used by       * * * */
/* tabipb() to interface with the treecode                   * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
int TreecodeInitialization(TABIPBparm *parm, struct Particles *particles)
{
    /* set up variables used in treecode */
    /* local variables*/
    int level;

    /* variables needed for reorder */
    double *temp_area, *temp_source;
    double *temp_normal[3];
    
    double xyz_limits[6];
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    int ierr;
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
    

/* Call SETUP to allocate arrays for Taylor expansions */
/* and setup global variables. Also, copy variables into global copy arrays. */
    s_Setup(xyz_limits, particles);
    
    int numnodes, numleaves, max_depth;
    
    TreeLinkedList_Construct(&s_tree_root, NULL, particles, 0, particles->num-1,
                s_max_per_leaf, xyz_limits, &numnodes, &numleaves,
                &s_min_level, &s_max_level, &max_depth, 0);
    
    if (rank == 0) {
        printf("Created tree for %d particles with max %d per node.\n\n",
               s_numpars, s_max_per_leaf);
    }
    
    reorder(particles->order.begin(), particles->order.end(), particles->nx.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->ny.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->nz.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->area.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->source_term.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->source_term.begin() + s_numpars);
    
    s_particle_position_x = particles->x.data();
    s_particle_position_y = particles->y.data();
    s_particle_position_z = particles->z.data();
    s_particle_normal_x = particles->nx.data();
    s_particle_normal_y = particles->ny.data();
    s_particle_normal_z = particles->nz.data();
    s_particle_area = particles->area.data();
    s_source_term = particles->source_term.data();
    
    s_target_charge_   = (double *)malloc(s_numpars * sizeof(double));
    s_target_charge_dx = (double *)malloc(s_numpars * sizeof(double));
    s_target_charge_dy = (double *)malloc(s_numpars * sizeof(double));
    s_target_charge_dz = (double *)malloc(s_numpars * sizeof(double));
    s_source_charge_   = (double *)malloc(s_numpars * sizeof(double));
    s_source_charge_dx = (double *)malloc(s_numpars * sizeof(double));
    s_source_charge_dy = (double *)malloc(s_numpars * sizeof(double));
    s_source_charge_dz = (double *)malloc(s_numpars * sizeof(double));

    return 0;
}
/**********************************************************/


/********************************************************/
int TreecodeFinalization(struct Particles *particles)
{
    double *temp_area, *temp_source, *temp_xvct;
    double *temp_normal[3], *temp_position[3];
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    int ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

/***********reorder particles*************/
    
    reorder(particles->order.begin(), particles->order.end(), particles->nx.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->ny.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->nz.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->area.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->source_term.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->source_term.begin() + particles->num);
    reorder(particles->order.begin(), particles->order.end(), particles->xvct.begin());
    reorder(particles->order.begin(), particles->order.end(), particles->xvct.begin() + particles->num);

    
/***********treecode_initialization*******/

    free(s_target_charge_);
    free(s_target_charge_dx);
    free(s_target_charge_dy);
    free(s_target_charge_dz);
    free(s_source_charge_);
    free(s_source_charge_dx);
    free(s_source_charge_dy);
    free(s_source_charge_dz);
    
/***********clean tree structure**********/

    TreeLinkedList_Free(&s_tree_root);

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
    double temp_x, temp_area;
    double temp_charge[4];
    double pre1, pre2;
    double peng[2], peng_old[2];
    double *tpoten_temp;
    
    int particles_per_process;
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    int ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif
    
    //make_vector(tpoten_temp, 2 * s_numpars);
    tpoten_temp = (double *)malloc(s_numpars * 2 * sizeof(double));
    memcpy(tpoten_temp, tpoten, 2 * s_numpars * sizeof(double));
    memset(tpoten, 0, 2 * s_numpars * sizeof(double));
    
    s_ComputePBKernel(tpoten_old);
    
  /* Generate the moments if not allocated yet */
    s_ComputeAllMoments(s_tree_root, 1);

    pre1 = 0.50 * (1.0 + s_eps);
    pre2 = 0.50 * (1.0 + 1.0/s_eps);
    
    particles_per_process = s_numpars / num_procs;

    for (int ii = 0; ii <= particles_per_process; ii++) {
        int i = ii * num_procs + rank;
        
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
        
            temp_charge[0] = s_target_charge_  [i];
            temp_charge[1] = s_target_charge_dx[i];
            temp_charge[2] = s_target_charge_dy[i];
            temp_charge[3] = s_target_charge_dz[i];

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

    free(tpoten_temp);

    s_RemoveMoments(s_tree_root);

    return 0;
}
/**********************************************************/


/**********************************************************/
int psolve(double *z, double *r)
{
/* r as original while z as scaled */

        double scale1 = 0.5 * (1.0 + s_eps);
        double scale2 = 0.5 * (1.0 + 1.0/s_eps);

        for (int i = 0; i < s_numpars; i++) {
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

    int idx = 0, nrow = 0, nrow2, ibeg = 0, iend = 0;
    int *ipiv;
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
    columnMajorA = (double *)malloc(4*s_max_per_leaf*s_max_per_leaf*sizeof(double));
    rhs          = (double *)malloc(2*s_max_per_leaf*sizeof(double));
    ipiv         = (int *)malloc(2*s_max_per_leaf*sizeof(int));

    while (idx < s_numpars) {
        leaflength(s_tree_root, idx, &nrow);
        nrow2 = nrow*2;
        ibeg  = idx;
        iend  = idx + nrow - 1;

        memset(columnMajorA, 0, nrow2*nrow2*sizeof(double));
        memset(ipiv, 0, nrow2*sizeof(int));
        memset(rhs, 0, nrow2*sizeof(double));

        for (int i = ibeg; i <= iend; i++) {
            tp[0] = s_particle_position_x[i];
            tp[1] = s_particle_position_y[i];
            tp[2] = s_particle_position_z[i];
            tq[0] = s_particle_normal_x[i];
            tq[1] = s_particle_normal_y[i];
            tq[2] = s_particle_normal_z[i];

            for (int j = ibeg; j < i; j++) {
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

                columnMajorA[(j-ibeg)*nrow2 + i-ibeg] = -L1*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i-ibeg] = -L2*area;
                columnMajorA[(j-ibeg)*nrow2 + i+nrow-ibeg] = -L3*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i+nrow-ibeg] = -L4*area;
            }

            columnMajorA[(i-ibeg)*nrow2 + i-ibeg] = pre1;
            columnMajorA[(i+nrow-ibeg)*nrow2 + i+nrow-ibeg] = pre2;

            for (int j = i+1; j <= iend; j++) {
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

                columnMajorA[(j-ibeg)*nrow2 + i-ibeg] = -L1*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i-ibeg] = -L2*area;
                columnMajorA[(j-ibeg)*nrow2 + i+nrow-ibeg] = -L3*area;
                columnMajorA[(j+nrow-ibeg)*nrow2 + i+nrow-ibeg] = -L4*area;
            }
        }

        for (int i = 0; i < nrow; i++) {
            rhs[i] = r[i+ibeg];
            rhs[i+nrow] = r[i+ibeg+s_numpars];
        }

        // Apple Accelerate implementation of LAPACK LU decomposition
        //nrhs = 1;
        //dgesv_(&nrow2, &nrhs, columnMajorA, &nrow2, ipiv, rhs, &nrow2, &info);

        // LAPACKE implementation of LAPACK LU decomposition
        LAPACKE_dgesv(LAPACK_COL_MAJOR, nrow2, 1, columnMajorA, nrow2, ipiv, rhs, nrow2);

        for (int i = 0; i < nrow; i++) {
            z[i+ibeg] = rhs[i];
            z[i+ibeg+s_numpars] = rhs[i+nrow];
        }

        idx += nrow;
    }

    //free_matrix(matrixA);
    free(columnMajorA);
    free(rhs);
    free(ipiv);
  
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
 
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    int ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif
    
    if (rank == 0) {
        printf("Setting up arrays for Taylor expansion...\n");
    }

    tt = (double *)malloc(s_torder_lim * sizeof(double));
    ww = (double *)malloc(s_torder_lim * sizeof(double));
    
    /* initializing array for Chev points */
    for (int i = 0; i < s_torder_lim; i++)
        tt[i] = cos(i * M_PI / s_order);
    
    ww[0] = 0.25 * (s_order*s_order/3.0 + 1.0/6.0);
    ww[s_order] = -ww[0];
    
    for (int i = 1; i < s_order; i++) {
        double xx = i * M_PI / s_order;
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
    }
    
    xyz_limits[0] = *std::min_element(particles->x.begin(), particles->x.end());
    xyz_limits[1] = *std::max_element(particles->x.begin(), particles->x.end());

    xyz_limits[2] = *std::min_element(particles->y.begin(), particles->y.end());
    xyz_limits[3] = *std::max_element(particles->y.begin(), particles->y.end());

    xyz_limits[4] = *std::min_element(particles->z.begin(), particles->z.end());
    xyz_limits[5] = *std::max_element(particles->z.begin(), particles->z.end());
    
    particles->order.resize(particles->num);
    std::iota(particles->order.begin(), particles->order.end(), 0);

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputePBKernel(double *phi)
{

    for (int i = 0; i < s_numpars; i++) {

        s_target_charge_  [i] = ONE_OVER_4PI;
        s_target_charge_dx[i] = ONE_OVER_4PI * s_particle_normal_x[i];
        s_target_charge_dy[i] = ONE_OVER_4PI * s_particle_normal_y[i];
        s_target_charge_dz[i] = ONE_OVER_4PI * s_particle_normal_z[i];
        
        s_source_charge_  [i] = s_particle_area[i] * phi[s_numpars+i];
        s_source_charge_dx[i] = s_particle_normal_x[i] * s_particle_area[i] * phi[i];
        s_source_charge_dy[i] = s_particle_normal_y[i] * s_particle_area[i] * phi[i];
        s_source_charge_dz[i] = s_particle_normal_z[i] * s_particle_area[i] * phi[i];

    }

    return 0;
}
/********************************************************/


/********************************************************/
static int s_ComputeAllMoments(struct TreeLinkedListNode *p, int ifirst)
{
/* REMOVE_NODE recursively removes each node from the tree and deallocates
 * its memory for MS array if it exits. */

    if (p->exist_ms == 0 && ifirst == 0) {
        //make_matrix(p->ms, 4, s_torder3);
        p->ms[0] = (double *)malloc(s_torder3 * sizeof(double));  
        p->ms[1] = (double *)malloc(s_torder3 * sizeof(double)); 
        p->ms[2] = (double *)malloc(s_torder3 * sizeof(double)); 
        p->ms[3] = (double *)malloc(s_torder3 * sizeof(double));  
        p->tx = (double *)malloc(s_torder_lim * sizeof(double));
        p->ty = (double *)malloc(s_torder_lim * sizeof(double));
        p->tz = (double *)malloc(s_torder_lim * sizeof(double));
        s_ComputeMoments(p);
        p->exist_ms = 1;
    }

    if (p->num_children > 0) {
        for (int i = 0; i < p->num_children; i++) {
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

    double dj[s_torder_lim];
    double a1i[s_torder_lim], a2j[s_torder_lim], a3k[s_torder_lim];
    double w1i[s_torder_lim];
    double summ[4][s_torder3];
    
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < s_torder3; j++) {
            p->ms[i][j] = 0.0;
        }
    }
    
    double* xibeg = &(s_particle_position_x[p->ibeg]);
    double* yibeg = &(s_particle_position_y[p->ibeg]);
    double* zibeg = &(s_particle_position_z[p->ibeg]);
    
    double x0 = p->x_min;
    double x1 = p->x_max;
    double y0 = p->y_min;
    double y1 = p->y_max;
    double z0 = p->z_min;
    double z1 = p->z_max;
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < s_torder3; j++) {
            summ[i][j] = 0.0;
        }
    }
    
    for (int i = 0; i < s_torder_lim; i++) {
        p->tx[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        p->ty[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        p->tz[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    
    dj[0] = 0.5;
    dj[s_order] = 0.5;
    for (int j = 1; j < s_order; j++)
        dj[j] = 1.0;
    
    for (int j = 0; j < s_torder_lim; j++)
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];

    for (int i = 0; i < p->numpar; i++) {
    
        double sumA1 = 0.0;
        double sumA2 = 0.0;
        double sumA3 = 0.0;
        
        int a1exactind = -1;
        int a2exactind = -1;
        int a3exactind = -1;
    
        double xx = xibeg[i];
        double yy = yibeg[i];
        double zz = zibeg[i];
        double qq_   = s_source_charge_  [p->ibeg+i];
        double qq_dx = s_source_charge_dx[p->ibeg+i];
        double qq_dy = s_source_charge_dy[p->ibeg+i];
        double qq_dz = s_source_charge_dz[p->ibeg+i];
        
        for (int j = 0; j < s_torder_lim; j++) {
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
            for (int j = 0; j < s_torder_lim; j++) a1i[j] = 0.0;
            a1i[a1exactind] = 1.0;
        }
        
        if (a2exactind > -1) {
            sumA2 = 1.0;
            for (int j = 0; j < s_torder_lim; j++) a2j[j] = 0.0;
            a2j[a2exactind] = 1.0;
        }
        
        if (a3exactind > -1) {
            sumA3 = 1.0;
            for (int j = 0; j < s_torder_lim; j++) a3k[j] = 0.0;
            a3k[a3exactind] = 1.0;
        }
        
        double Dd = 1.0 / (sumA1 * sumA2 * sumA3);

        int kk = -1;
        for (int k1 = 0; k1 < s_torder_lim; k1++) {
            for (int k2 = 0; k2 < s_torder_lim; k2++) {
                for (int k3 = 0; k3 < s_torder_lim; k3++) {
                    kk++;
                
                    double mom1 = a1i[k1] * a2j[k2] * a3k[k3] * Dd;
                    
                    summ[0][kk] += mom1 * qq_;
                    summ[1][kk] += mom1 * qq_dx;
                    summ[2][kk] += mom1 * qq_dy;
                    summ[3][kk] += mom1 * qq_dz;
                }
            }
        }
    }
    
    for (int j = 0; j < 4; j++)
        memcpy(p->ms[j], summ[j], s_torder3*sizeof(double));

    return 0;
}
/********************************************************/


/********************************************************/
static int s_RunTreecode(struct TreeLinkedListNode *p, double *tpoten_old, double tempq[4],
                         double peng[2])
{
  /* RunTreecode() is self recurrence function */
  
    double pengchild[2];


  /* determine DISTSQ for MAC test */
    double tx = p->x_mid - s_target_position[0];
    double ty = p->y_mid - s_target_position[1];
    double tz = p->z_mid - s_target_position[2];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);

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
            for (int i = 0; i < p->num_children; i++) {
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
static int s_ComputeTreePB(struct TreeLinkedListNode *p, double tempq[4], double peng[2])
{
    
    double pt_comp_ = 0.;
    double pt_comp_dx = 0.;
    double pt_comp_dy = 0.;
    double pt_comp_dz = 0.;
    
    double* __restrict__ cluster_x = p->tx;
    double* __restrict__ cluster_y = p->ty;
    double* __restrict__ cluster_z = p->tz;

    double* __restrict__ cluster_q_   = p->ms[0];
    double* __restrict__ cluster_q_dx = p->ms[1];
    double* __restrict__ cluster_q_dy = p->ms[2];
    double* __restrict__ cluster_q_dz = p->ms[3];
    
    double target_x = s_target_position[0];
    double target_y = s_target_position[1];
    double target_z = s_target_position[2];
    
    int ii = 0;
    for (int i = 0; i < s_torder_lim; i++) {
        for (int j = 0; j < s_torder_lim; j++) {
            for (int k = 0; k < s_torder_lim; k++) {

                double dx = target_x - cluster_x[i];
                double dy = target_y - cluster_y[j];
                double dz = target_z - cluster_z[k];

                double r2    = dx*dx + dy*dy + dz*dz;
                double r     = sqrt(r2);
                double rinv  = 1.0 / r;
                double r3inv = rinv  * rinv * rinv;
                double r5inv = r3inv * rinv * rinv;

                double expkr   =  exp(-s_kappa * r);
                double d1term  =  r3inv * expkr * (1. + (s_kappa * r)); 
                double d1term1 = -r3inv + d1term * s_eps;
                double d1term2 = -r3inv + d1term / s_eps;
                double d2term  =  r5inv * (-3. + expkr * (3. + (3. * s_kappa * r)
                                                       + (s_kappa * s_kappa * r2)));
                double d3term  =  r3inv * ( 1. - expkr * (1. + s_kappa * r));

                pt_comp_    += (rinv * (1. - expkr) * (cluster_q_  [ii])
                                          + d1term1 * (cluster_q_dx[ii] * dx
                                                     + cluster_q_dy[ii] * dy
                                                     + cluster_q_dz[ii] * dz));
                                        
                pt_comp_dx  += (cluster_q_  [ii]  * (d1term2 * dx)
                             - (cluster_q_dx[ii]  * (dx * dx * d2term + d3term)
                             +  cluster_q_dy[ii]  * (dx * dy * d2term)
                             +  cluster_q_dz[ii]  * (dx * dz * d2term)));
                             
                pt_comp_dy  += (cluster_q_  [ii]  *  d1term2 * dy
                             - (cluster_q_dx[ii]  * (dx * dy * d2term)
                             +  cluster_q_dy[ii]  * (dy * dy * d2term + d3term)
                             +  cluster_q_dz[ii]  * (dy * dz * d2term)));
                             
                pt_comp_dz  += (cluster_q_  [ii]  *  d1term2 * dz
                             - (cluster_q_dx[ii]  * (dx * dz * d2term)
                             +  cluster_q_dy[ii]  * (dy * dz * d2term)
                             +  cluster_q_dz[ii]  * (dz * dz * d2term + d3term)));

                ii++;
            }
        }
    }
            
    peng[0] = tempq[0] * pt_comp_;
    peng[1] = tempq[1] * pt_comp_dx  +  tempq[2] * pt_comp_dy  +  tempq[3] * pt_comp_dz;

    return 0;
}
/********************************************************/



/********************************************************/
static int s_ComputeDirectPB(int ibeg, int iend,
                             double *tpoten_old, double peng[2])
{
  /* COMPF_DIRECT directly computes the force on the current target
 * particle determined by the global variable s_target_position.*/
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
    
    for (int j = ibeg; j < iend+1; j++) {
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
 
    if (p->exist_ms == 1) {
        free(p->ms[0]);
        free(p->ms[1]);
        free(p->ms[2]);
        free(p->ms[3]);
        free(p->tx);
        free(p->ty);
        free(p->tz);
        p->exist_ms = 0;
    }

    if (p->num_children > 0) {
        for (int i = 0; i < p->num_children; i++)
            s_RemoveMoments(p->child[i]);
    }

    return 0;
}
/********************************************************/


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* PRIVATE functions for preconditioning                     * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
static void leaflength(struct TreeLinkedListNode *p, int idx, int *nrow)
{
/* find the leaf length */

    if (idx == p->ibeg && p->num_children == 0) {
        *nrow = p->numpar;
    } else {
        if (p->num_children != 0) {
            for (int i = 0; i < p->num_children; i++)
                leaflength(p->child[i], idx, nrow);
        }
     }

}
/**********************************************************/
