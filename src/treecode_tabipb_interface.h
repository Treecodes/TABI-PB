/*
 * C header file for treecode interface with tabipb
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

#ifndef H_TREECODE_TABIPB_INTERFACE_H
#define H_TREECODE_TABIPB_INTERFACE_H

#include "TABIPBstruct.h"
#include "particle_struct.h"


/* functions used by tabipb() to interface with treecode */
int TreecodeInitialization(TABIPBparm *parm, int nface,
                           TreeParticles *particles);

int TreecodeFinalization(TreeParticles *particles);

#endif /* H_TREECODE_TABIPB_INTERFACE_H */
