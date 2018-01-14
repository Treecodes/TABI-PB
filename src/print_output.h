/*
 * C header for routines for printing tabipb output
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

#include "TABIPBstruct.h"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* output printing functions                                 * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int OutputPrint(TABIPBvars *vars);

int OutputVTK(TABIPBparm *parm, TABIPBvars *vars);
