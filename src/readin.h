/**************************************************************************
* FILE NAME: readin.h                                                     *
*                                                                         *
* PURPOSE: header for readin routine responsible for meshing surface and  *
*          and reading in surface data, called by tabipb routine          *
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
* 01/12/2018  Leighton Wilson   Created                                   *
*                                                                         *
**************************************************************************/

#ifndef H_READIN_H
#define H_READIN_H

#include "TABIPBstruct.h"
#include "particle_struct.h"

/* function to read in molecule information */
int Readin(TABIPBparm *parm, TABIPBvars *vars);

#endif /* H_READIN_H */
