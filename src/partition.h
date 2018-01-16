/*
 * C header for partition routine used by treecode in tabipb
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Created by Leighton Wilson, 01/15/2018
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* helper functions                                          * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef H_PARTITION_H
#define H_PARTITION_H

int Partition(double *a, double *b, double *c, int *indarr,
              int ibeg, int iend, double val);

#endif /* H_PARTITION_H */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
