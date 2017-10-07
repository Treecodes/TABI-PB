/*
 * C routine to interface tabipb with apbs
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Additional modifications and updates by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Works for Sphinx by Jiahui at 7/14/2016
 * Rebuild the architecture of wrapper by Jiahui at 6/30/2016
 * Build matrix free and nanoshaper by Leighton at 6/23/2016
 *
 */

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gl_variables.h"
#include "treecode.h"
#include "array.h"
#include "TABIPBstruct.h"

int tabipb(TABIPBparm *parm, TABIPBvars *vars) {
  /* Assemble the TABIPBparm out side this subroutine, and pass the three arryas */
  /* TABIPBparm a structure of parameters: file path, file name, density,
     probe radius, epsp, epsw, bulk_strength, treecode order, treecode maxparnode,
     treecode theta, mesh flag, and number of lines in pqr file */
  /* t_chrpos, t_atmchr and t_atmrad are three arrays of position, charges and
     raduis */

  /* variables local to main */
    int i, j, k;
    double s[3];
    double pot = 0.0; 
    double sum = 0.0; 
    double pot_temp = 0.0;
    double ptl; 
    double soleng = 0.0; 
    double couleng = 0.0;
    double t1, t2;

    extern void readin();
    extern double potential_molecule(double s[3]);
    extern int comp_source(TABIPBparm *parm, TABIPBvars *vars);
    extern int output_potential();

    /* variables used to compute potential solution */
    double units_para;
    double *chrptl;
    extern int comp_pot(TABIPBvars *vars, double *chrptl);

    /* variables used in treecode */
    extern int treecode_initialization(int order, int maxparnode, double theta);
    extern int treecode_finalization();

    /* GMRES related variables */
    static long int info;
    long int RESTRT, ldw, ldh, iter, N;
    double resid;

    extern int matvec(), psolve();
    extern int gmres_(long int *n, double *b, double *x, long int *restrt,
                      double *work, long int *ldw, double *h, long int *ldh,
                      long int *iter, double *resid, int (*matvec)(),
                      int (*psolve)(), long int *info);

    printf("\n Treecode order: %d", parm->order);
    printf("\n Max # of nodes: %d", parm->maxparnode);
    printf("\n    MAC (theta): %f\n", parm->theta);
    printf("\n      Mesh flag: %d\n", parm->mesh_flag);

    printf("\nSetting up the TABI input...\n");

  /***************constant*****************/
    pi = 3.14159265358979324;
    one_over_4pi = 0.079577471545948;
    kcal2j = 4.184;
    bulk_coef = 2529.12179861515279; /* constant without temperature */
    units_coef = 332.0716 * kcal2j;
    eps = parm->epsw/parm->epsp;
    kappa2 = bulk_coef * parm->bulk_strength / parm->epsw / parm->temp;
    kappa = sqrt(kappa2);

  /*read charge coodinates, charges and radius*/
  //chrpos = &vars->chrpos[0];
  //atmchr = &vars->atmchr[0];
  //atmrad = &vars->atmrad[0];

    readin(parm, vars);

    bvct = (double *) calloc(2*nface, sizeof(double));
    comp_source(parm, vars);

    /* set up treecode */
    treecode_initialization(parm->order, parm->maxparnode, parm->theta);

    /* parameters for GMRES */
    RESTRT = 10;
    N = 2*nface;
    ldw = N;
    ldh = RESTRT + 1;
    iter = 100;
    resid = 1e-4;
    xvct = (double *) calloc(N, sizeof(double));

    work = (double *) calloc(ldw * (RESTRT + 4), sizeof(double));
    h = (double *) calloc(ldh * (RESTRT + 2), sizeof(double));

    gmres_(&N, bvct, xvct, &RESTRT, work, &ldw, h, &ldh, &iter,
           &resid, &matvec, psolve, &info);

    /* the solvation energy computation */
    units_para = 2.0 * units_coef * pi;

    chrptl = (double *) calloc(nface, sizeof(double));
    comp_pot(vars, chrptl);

    soleng = 0.0;
    couleng = 0.0;

    double r[3], diff[3], dist;

    for (i = 0; i < nface; i++)
        soleng += chrptl[i];

    for (i = 0; i < natm; i++) {
//  for (i = 0; i < 10; i++){
        r[0] = vars->chrpos[3*i];
        r[1] = vars->chrpos[3*i + 1];
        r[2] = vars->chrpos[3*i + 2];
        for (j = i+1; j < natm; j++){
            diff[0] = r[0] - vars->chrpos[3*j];
            diff[1] = r[1] - vars->chrpos[3*j + 1];
            diff[2] = r[2] - vars->chrpos[3*j + 2];
            dist = sqrt(diff[0]*diff[0]
                 + diff[1]*diff[1]
                 + diff[2]*diff[2]);
            couleng += 1 / parm->epsp / dist * vars->atmchr[i] * vars->atmchr[j];
        }
        //printf("the couleng is %f,%f\n",couleng,dist);
    }

    soleng = soleng * units_para;
    couleng = couleng * units_coef;
    vars->soleng = soleng;
    vars->couleng = couleng;

    output_potential(vars);

    free_matrix(vars->extr_v);
    free_matrix(vert);
    free_matrix(snrm);
    free_matrix(face);
    free_matrix(vars->extr_f);

    free(tr_xyz);
    free(tr_q);

    free(tr_area);
    free(bvct);
    free(xvct);

    /* dellocate treecode variables */
    treecode_finalization();

    return 0;

}

/************************************/
int psolve(double *z, double *r) {
/* r as original while z as scaled */

        int i;
        double scale1, scale2;
        scale1 = 0.5 * (1.0 + eps);
        scale2 = 0.5 * (1.0 + 1.0/eps);

        for (i = 0; i < nface; i++) {
                z[i] = r[i]/scale1;
                z[i + nface] = r[i + nface]/scale2;
        }

        return 0;
}

/************************************/
int comp_source(TABIPBparm *parm, TABIPBvars *vars) {
/* this computes the source term where
 * S1=sum(qk*G0)/e1 S2=sim(qk*G0')/e1 */

/* local variables */
    int i, j;
    double sumrs, cos_theta, irs, G0, G1, tp1;
    double r_s[3];

/* global variables located in gl_variables.h 
        double *tr_xyz;
        double *tr_q;  
        double *bvct;
 */

    for (i = 0; i < nface; i++) {
        bvct[i] = 0.0;
        bvct[i + nface] = 0.0;
        for (j = 0; j < natm; j++) {

  /* r_s = distance of charge position to triangular */
            r_s[0] = vars->chrpos[3*j] - tr_xyz[3*i];
            r_s[1] = vars->chrpos[3*j + 1] - tr_xyz[3*i + 1];
            r_s[2] = vars->chrpos[3*j + 2] - tr_xyz[3*i + 2];
            sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];

  /* cos_theta = <tr_q,r_s>/||r_s||_2 */
            cos_theta = tr_q[3*i] * r_s[0] + tr_q[3*i + 1] * r_s[1]
                      + tr_q[3*i + 2] * r_s[2];
            irs = 1 / sqrt(sumrs);
            cos_theta = cos_theta * irs;

  /* G0 = 1/(4pi*||r_s||_2) */
            G0 = one_over_4pi;
            G0 = G0 * irs;

  /* G1 = cos_theta*G0/||r_s||_2 */
            tp1 = G0 * irs;
            G1 = cos_theta * tp1;

  /* update bvct */
            bvct[i] += vars->atmchr[j] * G0 / parm->epsp;
            bvct[nface + i] += vars->atmchr[j] * G1 / parm->epsp;
        }
    }

    return 0;
}

/************************************/

/************************************/
int comp_pot(TABIPBvars *vars, double *chrptl) {
  /* local variables */
    int i, j;
    double sumrs, irs, rs, G0, Gk, kappa_rs, exp_kappa_rs;
    double cos_theta, G1, G2, L1, L2, tp1, tp2;
    double r[3], v[3], s[3], r_s[3];

  /* global variables from gl_variables.h 
        double *xvct;
        double *tr_area;
        double *tr_xyz;
        double *tr_q
  */

    for (j = 0; j < nface; j++) {
        chrptl[j] = 0.0;

  /* r[] = tr_xyz[] & v[] = tr_q[] */
        r[0] = tr_xyz[3*j];
        r[1] = tr_xyz[3*j + 1];
        r[2] = tr_xyz[3*j + 2];

        v[0] = tr_q[j*3];
        v[1] = tr_q[j*3 + 1];
        v[2] = tr_q[j*3 + 2];

        for (i = 0; i < natm; i++) {
  /* s = chrpos[] & r_s = r[]-s[] */
            s[0] = vars->chrpos[3*i];
            s[1] = vars->chrpos[3*i + 1];
            s[2] = vars->chrpos[3*i + 2];

            r_s[0] = r[0] - s[0];
            r_s[1] = r[1] - s[1];
            r_s[2] = r[2] - s[2];

            sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
            rs = sqrt(sumrs);
            irs = 1 / rs;

            G0 = one_over_4pi;
            G0 = G0 * irs;
            kappa_rs = kappa * rs;
            exp_kappa_rs = exp(-kappa_rs);
            Gk = exp_kappa_rs * G0;

            cos_theta = (v[0]*r_s[0] + v[1]*r_s[1] + v[2]*r_s[2]) * irs;

            tp1 = G0 * irs;
            tp2 = (1.0 + kappa_rs) * exp_kappa_rs;

            G1 = cos_theta * tp1;
            G2 = tp2 * G1;

            L1 = G1 - eps * G2;
            L2 = G0 - Gk;

            chrptl[j] += vars->atmchr[i]
                         * (L1*xvct[j] + L2*xvct[nface+j]) * tr_area[j];
        }
    }

    return 0;
}


/************************************/
int output_potential(TABIPBvars *vars) {

    int i, j, k, jerr, nface_vert;
    double tot_length, loc_length, aa[3], dot_aa, para_temp, phi_star;
    int **ind_vert;
    double *xtemp, *xyz_temp;

    extern double maxval(double *, int);
    extern double minval(double *, int);

    nface_vert = 15; /* one vertex could have been involved
                        in at most 11 triangles, 15 is safe */
    para_temp = units_coef * 4 * pi;

    if ((xtemp = (double *) calloc(2 * nface, sizeof(double)))  == NULL) {
        printf("Error in allocating xtemp!\n");
    }

    if ((ind_vert = (int **) calloc(nface_vert, sizeof(int *)))  == NULL) {
        printf("Error in allocating vars->xvct!\n");
    }

    for (i = 0; i < nface_vert; i++){
        if ((ind_vert[i] = (int *) calloc(nspt, sizeof(int)))  == NULL) {
             printf("Error in allocating vars->xvct!\n");
        }
    }

    if ((vars->vert_ptl = (double *) calloc(nspt * 2, sizeof(double)))  == NULL) {
        printf("Error in allocating vars->xvct!\n");
    }
    if ((xyz_temp = (double *) calloc(3 * nface, sizeof(double)))  == NULL) {
        printf("Error in allocating vars->xvct!\n");
    }

  /* put things back */
    for (i = 0; i < nface; i++) {
        xtemp[orderarr[i]] = xvct[i];
        xtemp[orderarr[i] + nface] = xvct[i + nface];
        xyz_temp[orderarr[i]*3] = tr_xyz[i*3];
        xyz_temp[orderarr[i]*3 + 1] = tr_xyz[i*3 + 1];
        xyz_temp[orderarr[i]*3 + 2] = tr_xyz[i*3 + 2];
    }

    for (i = 0; i < nface; i++) {
        xvct[i] = xtemp[i];
        xvct[i + nface] = xtemp[i + nface];
        tr_xyz[i*3] = xyz_temp[i*3];
        tr_xyz[i*3 + 1] = xyz_temp[i*3 + 1];
        tr_xyz[i*3 + 2] = xyz_temp[i*3 + 2];
    }

    for (i = 0; i < nface; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < nface_vert - 1; k++) {
                if (ind_vert[k][face[j][i] - 1] == 0.0) {
                    ind_vert[k][face[j][i] - 1] = i + 1;
                    ind_vert[nface_vert - 1][face[j][i] - 1] += 1;
                    break;
                }
            }
        }
    }

    for (i = 0; i < nspt; i++) {
        tot_length = 0.0;
        for (j = 0; j < ind_vert[nface_vert - 1][i]; j++) {
      /* distance between vertices and centroid */
            aa[0] = tr_xyz[3 * (ind_vert[j][i]-1)] - vert[0][i];
            aa[1] = tr_xyz[3 * (ind_vert[j][i]-1) + 1] - vert[1][i];
            aa[2] = tr_xyz[3 * (ind_vert[j][i]-1) + 2] - vert[2][i];
            dot_aa = aa[0]*aa[0] + aa[1]*aa[1] + aa[2]*aa[2];
            loc_length = sqrt(dot_aa);

            vars->vert_ptl[i] += 1.0/loc_length*xvct[ind_vert[j][i]-1];
            vars->vert_ptl[i+nspt] += 1.0/loc_length*xvct[ind_vert[j][i]+nface-1];
            tot_length += 1.0/loc_length;
        }
        vars->vert_ptl[i] /= tot_length;
        vars->vert_ptl[i + nspt] /= tot_length;
    }

    for (i = 0; i < 2 * nface; i++)
        xvct[i] *= para_temp;

    for (i = 0; i < nspt; i++) {
        vars->vert_ptl[i] *= para_temp;
        vars->vert_ptl[i + nspt] *= para_temp;
    }

    if ((vars->xvct = (double *) calloc(2 * nface, sizeof(double))) == NULL) {
        printf("Error in allocating vars->xvct!\n");
    }
    memcpy(vars->xvct, xvct, 2 * nface * sizeof(double));

    vars->max_xvct = maxval(xvct, nface);
    vars->min_xvct = minval(xvct, nface);
    vars->max_der_xvct = maxval(xvct + nface, nface);
    vars->min_der_xvct = minval(xvct + nface, nface);

    vars->max_vert_ptl = maxval(vars->vert_ptl, nspt);
    vars->min_vert_ptl = minval(vars->vert_ptl, nspt);
    vars->max_der_vert_ptl = maxval(vars->vert_ptl + nspt, nspt);
    vars->min_der_vert_ptl = minval(vars->vert_ptl + nspt, nspt);

    vars->nface = nface;
    vars->nspt = nspt;
    make_matrix(vars->vert, 3, nspt);
    make_matrix(vars->snrm, 3, nspt);
    make_matrix(vars->face, 3, nface);
    for (i = 0; i < 3; i++)
        memcpy(vars->vert[i], vert[i], nspt * sizeof(double));
    for (i = 0; i < 3; i++)
        memcpy(vars->snrm[i], snrm[i], nspt * sizeof(double));
    for (i = 0; i < 3; i++)
        memcpy(vars->face[i], face[i], nface * sizeof(int));

    free(xtemp);

    for (i = 0; i < nface_vert; i++)
        free(ind_vert[i]);

    free(ind_vert);
    free(xyz_temp);

    return 0;
}


/************************************/
int output_print(TABIPBvars *vars)
{
    int i;

    printf("\nSolvation energy = %f kJ/mol", vars->soleng);
    printf("\nFree energy = %f kJ/mol\n\n", vars->soleng+vars->couleng);
    printf("The max and min potential and normal derivatives on elements area:\n");
    printf("potential %f %f\n", vars->max_xvct, vars->min_xvct);
    printf("norm derv %f %f\n\n", vars->max_der_xvct,
                                  vars->min_der_xvct);
    printf("The max and min potential and normal derivatives on vertices area:\n");
    printf("potential %f %f\n", vars->max_vert_ptl, vars->min_vert_ptl);
    printf("norm derv %f %f\n\n", vars->max_der_vert_ptl,
                                  vars->min_der_vert_ptl);

    FILE *fp = fopen("surface_potential.dat", "w");
    fprintf(fp, "%d %d\n", vars->nspt, vars->nface);

    for (i = 0; i < nspt; i++)
        fprintf(fp, "%d %f %f %f %f %f %f %f %f\n", i,
                vars->vert[0][i], vars->vert[1][i], vars->vert[2][i],
                vars->snrm[0][i], vars->snrm[1][i], vars->snrm[2][i],
                vars->vert_ptl[i], vars->vert_ptl[i + nspt]);

    for (i = 0; i < nface; i++)
        fprintf(fp, "%d %d %d\n", vars->face[0][i], vars->face[1][i], vars->face[2][i]);

    fclose(fp);

    return 0;
}


/************************************/
int output_vtk(TABIPBparm *parm, TABIPBvars *vars)
{
    char i_char1[20], i_char2[20], i_char3[20], nspt_str[20],
         nface_str[20], nface4_str[20];
    int i;
    double cal2j=4.184;

    sprintf(nspt_str, "%d", vars->nspt);
    sprintf(nface_str, "%d", vars->nface);
    sprintf(nface4_str, "%d", vars->nface * 4);

    sprintf(i_char1, "mesh flag: %d", parm->mesh_flag);

    FILE *fp = fopen("surface_potential.vtk", "w");

    fprintf(fp, "# vtk DataFile Version 1.0\n");
    fprintf(fp, "mesh for protein %s, with %s\n", parm->fname, i_char1);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n\n");

    fprintf(fp, "POINTS %s double\n", nspt_str);
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f %f %f\n", vars->vert[0][i], vars->vert[1][i],
                                  vars->vert[2][i]);
    }

    fprintf(fp, "POLYGONS %s %s\n", nface_str, nface4_str);
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "3 %d %d %d\n", vars->face[0][i] - 1, vars->face[1][i] - 1,
                                    vars->face[2][i] - 1);
    }

    fprintf(fp, "\nPOINT_DATA %s\n", nspt_str);
    fprintf(fp, "SCALARS PotentialVert double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nspt; i++) {
         fprintf(fp, "%f\n", cal2j*vars->vert_ptl[i]);
    }

    fprintf(fp, "SCALARS NormalPotentialVert double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f\n", cal2j*vars->vert_ptl[nspt + i]);
    }

    //if we want induced surface charges, we can multiply vertnorm by (1/eps + 1)
    fprintf(fp, "\nNORMALS VertNorms double\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f %f %f\n", vars->snrm[0][i], vars->snrm[1][i],
                                  vars->snrm[2][i]);
    }

    fprintf(fp, "\nCELL_DATA %s\n", nface_str);
    fprintf(fp, "SCALARS PotentialFace double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "%f\n", cal2j*vars->xvct[i]);
    }

    //if we want induced surface charges, we can multiply vertnorm by (1/eps + 1)
    fprintf(fp, "SCALARS NormalPotentialFace double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "%f\n", cal2j*vars->xvct[nface + i]);
    }

    fclose(fp);

    return 0;
}


/************************************/
int *matvec_direct(double *alpha, double *x, double *beta, double *y)
{
    int i, j;
    double pre1, pre2;
    double area, rs, irs, sumrs;
    double G0, kappa_rs, exp_kappa_rs, Gk;
    double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
    double G10, G20, G1, G2, G3, G4;
    double L1, L2, L3, L4;
    double tp[3], tq[3], sp[3], sq[3], r_s[3];
    double peng[2], peng_old[2];

    pre1 = 0.50 * (1.0 + eps); /* eps=80.0 a constant */
    pre2 = 0.50 * (1.0 + 1.0/eps); /* fdivide */

    for (i = 0; i < nface; i++) {

        tp[0] = tr_xyz[3*i];
        tp[1] = tr_xyz[3*i + 1];
        tp[2] = tr_xyz[3*i + 2];

        tq[0] = tr_q[3*i];
        tq[1] = tr_q[3*i + 1];
        tq[2] = tr_q[3*i + 2];

        peng[0] = 0.0;
        peng[1] = 0.0;

        for (j = 0; j < nface; j++) {

            if (j != i) {

                sp[0] = tr_xyz[3*j];
                sp[1] = tr_xyz[3*j + 1];
                sp[2] = tr_xyz[3*j + 2];

                sq[0] = tr_q[3*j];
                sq[1] = tr_q[3*j + 1];
                sq[2] = tr_q[3*j + 2];

                r_s[0] = sp[0] - tp[0];
                r_s[1] = sp[1] - tp[1];
                r_s[2] = sp[2] - tp[2];

                sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
                rs = sqrt(sumrs);
                irs = 1 / rs;
                G0 = one_over_4pi;
                G0 = G0 * irs;
                kappa_rs = kappa * rs;
                exp_kappa_rs = exp(-kappa_rs);
                Gk = exp_kappa_rs * G0;

                cos_theta = (sq[0]*r_s[0] + sq[1]*r_s[1]
                                          + sq[2]*r_s[2]) * irs;
                cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1]
                                           + tq[2]*r_s[2]) * irs;

                tp1 = G0 * irs;
                tp2 = (1.0 + kappa_rs) * exp_kappa_rs;

                G10 = cos_theta0 * tp1;
                G20 = tp2 * G10;

                G1 = cos_theta * tp1;
                G2 = tp2 * G1;

                dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
                G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs * tp1;
                G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
                L1 = G1 - eps*G2;
                L2 = G0 - Gk;
                L3 = G4 - G3;
                L4 = G10 - G20 / eps; /* fdivide */
                /* x involve first */
                peng_old[0] = x[j];
                peng_old[1] = x[j + nface];
                area = tr_area[j];
                peng[0] = peng[0] + (L1*peng_old[0] + L2*peng_old[1]) * area;
                peng[1] = peng[1] + (L3*peng_old[0] + L4*peng_old[1]) * area;
            } 
        }

    /* update the y value */
        y[i] = y[i] * *beta + (pre1*x[i] - peng[0]) * *alpha;
        y[nface+i] = y[nface+i] * *beta + (pre2*x[nface+i] - peng[1]) * *alpha;

    }

    return 0;
}
