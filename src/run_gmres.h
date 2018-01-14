/*
 * C header for routine for setting up and calling GMRes
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/12/2018
 */

#ifndef H_RUN_GMRES_H
#define H_RUN_GMRES_H

int RunGMRES(int nface, double *source_term, double *xvct);

#endif /* H_RUN_GMRES_H */
