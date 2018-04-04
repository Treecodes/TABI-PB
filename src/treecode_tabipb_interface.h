/**************************************************************************
* FILE NAME: treecode_tabipb_interface.h                                  *
*                                                                         *
* PURPOSE: header for functions that interface treecode with the primary  *
*          TABI-PB routine for tree initalization and finalization        *
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

#ifndef H_TREECODE_TABIPB_INTERFACE_H
#define H_TREECODE_TABIPB_INTERFACE_H

#include "TABIPBstruct.h"
#include "particle_struct.h"

/* functions used by tabipb() to interface with treecode */
int TreecodeInitialization(TABIPBparm *parm, TreeParticles *particles);

int TreecodeFinalization(TreeParticles *particles);

#endif /* H_TREECODE_TABIPB_INTERFACE_H */
