/**************************************************************************
* FILE NAME: tabipb.c                                                     *
*                                                                         *
* PURPOSE: Setup and run the TABI-PB solver                               *
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

#ifdef MPI_ENABLED
    #include <mpi.h>
#endif

#include "tabipb.h"
#include "treecode_tabipb_interface.h"
#include "run_gmres.h"
#include "utilities.h"

#include "global_params.h"
#include "array.h"
#include "TABIPBstruct.h"
#include "struct_particles.h"

/* internal functions */
static int s_ConstructTreeParticles(TABIPBvars *vars,
                                    struct Particles *particles);
static int s_ComputeSourceTerm(TABIPBparm *parm, TABIPBvars *vars,
                               struct Particles *particles);
static int s_OutputPotential(TABIPBvars *vars, struct Particles *particles);

static double s_ComputeSolvationEnergy(TABIPBparm *parm, TABIPBvars *vars,
                                       struct Particles *particles);
static double s_ComputeCoulombEnergy(TABIPBparm *parm, TABIPBvars *vars);


/********************************************************/
int TABIPB(TABIPBparm *parm, TABIPBvars *vars) {

    double energy_solvation;
    double energy_coulomb;
    struct Particles *particles = malloc(sizeof *particles);
    
    int rank = 0, num_procs = 1, ierr;
    long int iter;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

    if (rank == 0) {
        printf("\nSetting up the TABI input...\n");
    }

    /* set up needed TABIPBparm constants */
    parm->eps = parm->epsw/parm->epsp;
    parm->kappa2 = BULK_COEFF * parm->bulk_strength / parm->epsw / parm->temp;
    parm->kappa = sqrt(parm->kappa2);
    
    /*set up the particles struct used by treecode */
    s_ConstructTreeParticles(vars, particles);
    s_ComputeSourceTerm(parm, vars, particles);

    /* set up treecode */
    TreecodeInitialization(parm, particles);

    /* call GMRES */
    RunGMRES(particles->num, particles->source_term,
             parm->precond, particles->xvct, &iter);

    /* compute solvation and coulombic energy */
    energy_solvation = s_ComputeSolvationEnergy(parm, vars, particles);
    energy_coulomb = s_ComputeCoulombEnergy(parm, vars);

    vars->soleng = energy_solvation * UNITS_PARA;
    vars->couleng = energy_coulomb * UNITS_COEFF;
    vars->gmres_iter = (int)iter;

    /* deallocate treecode variables and reorder particles */
    TreecodeFinalization(particles);

    /* scale potential to proper units and copy to TABIPBvars */
    s_OutputPotential(vars, particles);

    /* deconstruct particles */
    free_vector(particles->x);
    free_vector(particles->y);
    free_vector(particles->z);
    
    free_vector(particles->nx);
    free_vector(particles->ny);
    free_vector(particles->nz);
    
    free_vector(particles->area);
    free_vector(particles->source_term);
    free_vector(particles->xvct);
    free(particles);
    
    return 0;
}
/********************************************************/



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* PRIVATE internal TABIPB functions                         * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**********************************************************/
static int s_ConstructTreeParticles(TABIPBvars *vars, struct Particles *particles)
{
/* function to construct particles used by tree */

    int i, j, k, ierr, idx[3];
    double r[3][3], sum = 0.0;
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif


    //NODE PATCH METHOD
    particles->num = vars->nspt;
    
    make_vector(particles->x, particles->num);
    make_vector(particles->y, particles->num);
    make_vector(particles->z, particles->num);
    
    make_vector(particles->nx, particles->num);
    make_vector(particles->ny, particles->num);
    make_vector(particles->nz, particles->num);
    
    make_vector(particles->area, particles->num);
    make_vector(particles->source_term, 2 * particles->num);
    make_vector(particles->xvct, 2 * particles->num);

    for (i = 0; i < vars->nspt; i++) {
        particles->x[i] = vars->vert[0][i];
        particles->y[i] = vars->vert[1][i];
        particles->z[i] = vars->vert[2][i];
        
        particles->nx[i] = vars->snrm[0][i];
        particles->ny[i] = vars->snrm[1][i];
        particles->nz[i] = vars->snrm[2][i];

        particles->area[i] = 0.0;
    }

    for (i = 0; i < vars->nface; i++) {
        for (j = 0; j < 3; j++) {
            idx[j] = vars->face[j][i];
        }

        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                r[j][k] = vars->vert[j][idx[k]-1];
            }
        }

        for (j = 0; j < 3; j++) {
            particles->area[idx[j]-1] += TriangleArea(r);
        }
    }



    for (i = 0; i < particles->num; i++) {
        particles->area[i] /= 3.0;
        sum += particles->area[i];
    }

    if (rank == 0) {
        printf("Total suface area = %.17f\n", sum);
        vars->surface_area = sum;
    }

    return 0;
}
/**********************************************************/


/**********************************************************/
static int s_ComputeSourceTerm(TABIPBparm *parm, TABIPBvars *vars,
                               struct Particles *particles)
{
/* this computes the source term where
 * S1=sum(qk*G0)/e1 S2=sim(qk*G0')/e1 */

/* local variables */
    int i, j, ii, faces_per_process, ierr;
    double sumrs, cos_theta, irs, G0, G1;
    double r_s[3];

    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

    faces_per_process = particles->num / num_procs;

    for (ii = 0; ii <= faces_per_process; ii++) {
        i = ii * num_procs + rank;
        
    if (i < particles->num) {
        particles->source_term[i] = 0.0;
        particles->source_term[i + particles->num] = 0.0;

        for (j = 0; j < vars->natm; j++) {

  /* r_s = distance of charge position to triangular */
            r_s[0] = vars->chrpos[3*j]     - particles->x[i];
            r_s[1] = vars->chrpos[3*j + 1] - particles->y[i];
            r_s[2] = vars->chrpos[3*j + 2] - particles->z[i];
            sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];

  /* cos_theta = <tr_q,r_s>/||r_s||_2 */
            cos_theta = particles->nx[i] * r_s[0]
                      + particles->ny[i] * r_s[1]
                      + particles->nz[i] * r_s[2];
            irs = 1 / sqrt(sumrs);
            cos_theta = cos_theta * irs;

  /* G0 = 1/(4pi*||r_s||_2) */
            G0 = ONE_OVER_4PI * irs;

  /* G1 = cos_theta*G0/||r_s||_2 */
            G1 = cos_theta * G0 * irs;

  /* update source term */
            particles->source_term[i] += vars->atmchr[j] * G0 / parm->epsp;
            particles->source_term[particles->num + i] += vars->atmchr[j] * G1 / parm->epsp;
        }
    }
    }

#ifdef MPI_ENABLED
    ierr = MPI_Allreduce(MPI_IN_PLACE, particles->source_term, 2 * particles->num,
                         MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    return 0;
}
/********************************************************/


/********************************************************/
static double s_ComputeSolvationEnergy(TABIPBparm *parm, TABIPBvars *vars,
                                       struct Particles *particles)
{
  /* local variables */
    int i, j, ii, ierr, atms_per_process;
    double sumrs, irs, rs, G0, Gk, kappa_rs, exp_kappa_rs;
    double cos_theta, G1, G2, L1, L2;
    double r[3], v[3], s[3], r_s[3];
    double *chrptl;
    double energy_solvation = 0.0;

    int rank = 0, num_procs = 1;

#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

    make_vector(chrptl, particles->num);
    atms_per_process = vars->natm / num_procs;
        
    for (j = 0; j < particles->num; j++) {

        chrptl[j] = 0.0;

        r[0] = particles->x[j];
        r[1] = particles->y[j];
        r[2] = particles->z[j];

        v[0] = particles->nx[j];
        v[1] = particles->ny[j];
        v[2] = particles->nz[j];


        for (ii = 0; ii <= atms_per_process; ii++) {
            i = ii * num_procs + rank;

        if (i < vars->natm) {

            s[0] = vars->chrpos[3*i];
            s[1] = vars->chrpos[3*i + 1];
            s[2] = vars->chrpos[3*i + 2];

            r_s[0] = r[0] - s[0];
            r_s[1] = r[1] - s[1];
            r_s[2] = r[2] - s[2];

            sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
            rs = sqrt(sumrs);
            irs = 1 / rs;

            G0 = ONE_OVER_4PI * irs;
            kappa_rs = parm->kappa * rs;
            exp_kappa_rs = exp(-kappa_rs);
            Gk = exp_kappa_rs * G0;

            cos_theta = (v[0]*r_s[0] + v[1]*r_s[1] + v[2]*r_s[2]) * irs;

            G1 = G0 * cos_theta * irs;
            G2 = G1 * (1.0 + kappa_rs) * exp_kappa_rs;

            L1 = G1 - parm->eps * G2;
            L2 = G0 - Gk;

            chrptl[j] += vars->atmchr[i]
                       * (L1*particles->xvct[j]
                       + L2*particles->xvct[particles->num + j])
                       * particles->area[j];
        }
        }
    }
    
#ifdef MPI_ENABLED
    ierr = MPI_Allreduce(MPI_IN_PLACE, chrptl, particles->num,
                         MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    for (i = 0; i < particles->num; i++) {
        energy_solvation += chrptl[i];
    }
    
    free_vector(chrptl);

    return (energy_solvation);
}
/********************************************************/


/********************************************************/
static double s_ComputeCoulombEnergy(TABIPBparm *parm, TABIPBvars *vars)
{
    int i, j;
    double r[3], diff[3], dist;
    double energy_coulomb = 0.0;

    for (i = 0; i < vars->natm; i++) {
        r[0] = vars->chrpos[3*i];
        r[1] = vars->chrpos[3*i + 1];
        r[2] = vars->chrpos[3*i + 2];
        
        for (j = i+1; j < vars->natm; j++) {
            diff[0] = r[0] - vars->chrpos[3*j];
            diff[1] = r[1] - vars->chrpos[3*j + 1];
            diff[2] = r[2] - vars->chrpos[3*j + 2];
            dist = sqrt(diff[0]*diff[0]
                 + diff[1]*diff[1]
                 + diff[2]*diff[2]);
            energy_coulomb += 1 / parm->epsp / dist
                            * vars->atmchr[i] * vars->atmchr[j];
        }
    }

    return (energy_coulomb);
}
/********************************************************/


/********************************************************/
static int s_OutputPotential(TABIPBvars *vars, struct Particles *particles)
{
    int i, j;
    double para_temp;

    para_temp = UNITS_COEFF * 4 * PI;
    
    memcpy(vars->vert_ptl, particles->xvct, 2 * vars->nspt * sizeof(double));

    for (i = 0; i < vars->nface; i++) {
        vars->xvct[i] = 0.0;
        vars->xvct[vars->nface + i] = 0.0;
        for (j = 0; j < 3; j++) {
            vars->xvct[i] += vars->vert_ptl[vars->face[j][i]-1];
            vars->xvct[vars->nface + i] += vars->vert_ptl[vars->nspt + vars->face[j][i]-1];
        }
        vars->xvct[i] /= 3.0;
        vars->xvct[vars->nface + i] /= 3.0;
    }

    for (i = 0; i < 2 * vars->nface; i++)
        vars->xvct[i] *= para_temp;

    for (i = 0; i < 2 * vars->nspt; i++)
        vars->vert_ptl[i] *= para_temp;

    vars->max_xvct = MaxVal(vars->xvct, vars->nface);
    vars->min_xvct = MinVal(vars->xvct, vars->nface);
    vars->max_der_xvct = MaxVal(vars->xvct + vars->nface, vars->nface);
    vars->min_der_xvct = MinVal(vars->xvct + vars->nface, vars->nface);

    vars->max_vert_ptl = MaxVal(vars->vert_ptl, vars->nspt);
    vars->min_vert_ptl = MinVal(vars->vert_ptl, vars->nspt);
    vars->max_der_vert_ptl = MaxVal(vars->vert_ptl + vars->nspt, vars->nspt);
    vars->min_der_vert_ptl = MinVal(vars->vert_ptl + vars->nspt, vars->nspt);

    return 0;
}
/********************************************************/
