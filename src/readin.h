/*
 * C header to routine to mesh surfaces and read in surface data for tabipb
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


#ifndef H_READIN_H
#define H_READIN_H

#include "TABIPBstruct.h"
#include "particle_struct.h"

/* function to read in molecule information */
int Readin(TABIPBparm *parm, TABIPBvars *vars);

#endif /* H_READIN_H */
