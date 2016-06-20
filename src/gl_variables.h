/*
 * C header file for global variables of tabipb
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
 * Last modified by Leighton Wilson, 06/20/2016
 */

#ifndef H_GLVARS_H
#define H_GLVARS_H


/*constant variables */
double pi, one_over_4pi, bulk_coef, units_coef, epsw, epsp, eps;
double bulk_strength, kappa2, kappa;


/*global scalar variables*/
int nface, nspt, natm, nchr;


/*dynamic allocated variables*/
int **extr_v; //[3][nspt]
int **extr_f; //[2][nface]
int **face, **face_copy;//[3][nface]

double **vert, **snrm; //[3][nspt];
double *tr_xyz, *tr_q; //[3][nface]
double *tr_area, *bvct, *xvct; //[nface];
double **atmpos; //[3][natm/nchr];
double *atmrad, *atmchr, *chrpos; //[natm/nchr];

double *work, *h;


/*device pointers*/
double *h_pot;
double *dev_xp, *dev_yp, *dev_zp, *dev_q, *dev_pot;


#endif /* H_GLVARS_H */
