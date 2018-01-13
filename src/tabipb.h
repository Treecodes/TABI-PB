/*
 * C header file for tabipb routines
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/11/2018
 * Created by Leighton Wilson, 01/11/2018
 */

#ifndef H_TABIPB_H
#define H_TABIPB_H

#include "TABIPBstruct.h"


int tabipb(TABIPBparm *parm, TABIPB *vars);
    //Readin
    //ComputeSourceTerm
    //TreecodeInitialization
    //ComputePotential
    //OutputPotential
    //TreecodeFinalization

int ComputeSourceTerm(TABIPBparm *parm, TABIPBvars *vars);

int ComputePotential(TABIPBvars *vars, double *chrptl);

int OutputPotential(TABIPBvars *vars);
    //maxval
    //minval

int OutputPrint(TABIPBvars *vars);

int OutputVTK(TABIPBparm *parm, TABIPBvars *vars);

#endif /* H_TABIPB_H */
