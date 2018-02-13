/**************************************************************************
* FILE NAME: print_output.h                                               *
*                                                                         *
* PURPOSE: header for printing routines called by main (when running      *
*          standalone) or the APBS wrapper after returning from the       *
*          primary tabipb routine                                         *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 02/10/2018  Leighton Wilson   Separated out DAT creation file           *
* 01/12/2018  Leighton Wilson   Created                                   *
*                                                                         *
**************************************************************************/

#ifndef H_PRINT_OUTPUT_H
#define H_PRINT_OUTPUT_H

#include "TABIPBstruct.h"

int OutputPrint(TABIPBvars *vars);
int OutputDAT(TABIPBparm *parm, TABIPBvars *vars);
int OutputVTK(TABIPBparm *parm, TABIPBvars *vars);

#endif /* H_PRINT_OUTPUT_H */
