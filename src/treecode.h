/*
 * C header file for treecode routine of tabipb
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 06/20/2016
 */

#ifndef H_TREECODE_H
#define H_TREECODE_H

#include "TABIPBstruct.h"
#include "particle_struct.h"


/* functions used by tabipb() to interface with treecode */
int TreecodeInitialization(TABIPBparm *parm, int nface,
                           TreeParticles *particles);

int TreecodeFinalization(TreeParticles *particles);


/* functions used by GMRes */
int matvec(double *alpha, double *tpoten_old,
           double *beta, double *tpoten);

int psolve(double *z, double *r);


#endif /* H_TREECODE_H */
