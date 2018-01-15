/*
 * C header for helper routines for treecode in tabipb
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* helper functions                                          * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef H_UTILITY_ROUTINES_H
#define H_UTILITY_ROUTINES_H

double MinVal(double *variables, int number);

double MaxVal(double *variables, int number);

int Partition(double *a, double *b, double *c, int *indarr,
              int ibeg, int iend, double val);

#endif /* H_UTILITY_ROUTINES_H */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
